#include "discover.h"
#include "akhal/gfa.h"
#include "akhal/gaf.h"

// Read name -> breakpoint index. Keys are strdup'd (see parse_gaf) because the
// akhal GAF reader owns/reuses the record strings, so the names must be copied
// to remain valid as hash keys. They are freed in discover().
KHASHL_MAP_INIT(KH_LOCAL, mapstr_t, mapstr, const char*, uint32_t, kh_hash_str, kh_eq_str)

KSEQ_INIT(gzFile, gzread)

/**
 * @brief Walk a GAF path string, collecting the visited node ids.
 *
 * Validates that every node exists in the graph and that the read maps to a
 * single strand.
 *
 * @param readName Read name (for diagnostics).
 * @param readPath GAF path string (e.g. ">1>2<3").
 * @param nodes    Output: newly allocated array of node ids (caller frees).
 * @param g        The graph.
 * @return Number of nodes on success, or 0 on error.
 */
int parse_nodes(const char *readName, const char *readPath, uint64_t **nodes, const gfa_t *g) {

    uint64_t *temp_nodes = (uint64_t *)malloc(sizeof(uint64_t) * MAX_CIGAR);
    validate(temp_nodes)
    int node_size = 0, path_index, fwd_strand_count = 0, rev_strand_count = 0;
    const char *path = readPath; uint64_t id = 0; char strand = '>';

    while ((path_index = next_node(path, &id, &strand))) {
        path += path_index;
        if (gfa_idx(g, id) < 0) {
            fprintf(stderr, "[ERROR] Segment %lu in read %s not found.\n", id, readName);
            free(temp_nodes);
            return 0;
        }
        temp_nodes[node_size++] = id;

        if (strand == '>') fwd_strand_count++;
        else rev_strand_count++;
        // Validation: does read align consistently to one strand
        if (fwd_strand_count && rev_strand_count) {
            fprintf(stderr, "[ERROR] Read %s aligned in mixed strands.\n", readName);
            free(temp_nodes);
            return 0;
        }
    }

    if (node_size == 0) { free(temp_nodes); return 0; }  // don't leak on empty path
    *nodes = temp_nodes;
    return node_size;
}

/**
 * @brief Read a GAF file and derive per-read breakpoints against the graph.
 *
 * Uses the akhal GAF reader; a segment is treated as reference iff it belongs
 * to a path (ref_name != NULL). CIGAR handling stays svtrek-native (int-coded
 * parse_cigar).
 *
 * @param file_path       Path to the GAF file.
 * @param g               The graph.
 * @param gaf_breakpoints Output: newly allocated breakpoint array.
 * @param h2              Read name -> breakpoint index map (keys strdup'd here).
 * @return Number of breakpoints, or 0 if the file could not be opened.
 */
int parse_gaf(const char* file_path, const gfa_t *g, gaf_breakpoint **gaf_breakpoints, mapstr_t *h2) {

    gaf_reader_t *reader = gaf_open(file_path);
    if (!reader) { *gaf_breakpoints = NULL; return -1; }   // -1 distinguishes open failure

    int temp_breakpoint_size = 0, temp_breakpoint_capacity = 1000000;
    gaf_breakpoint *temp_breakpoints = (gaf_breakpoint *)malloc(temp_breakpoint_capacity * sizeof(gaf_breakpoint));
    validate(temp_breakpoints)

    gaf_rec_t rec;
    gaf_rec_init(&rec);

    while (gaf_read1(reader, &rec) == 1) {

        // Local, mutable copies of the fields we read/rewrite.
        const char *readName = rec.qname;
        int readLen   = (int)rec.qlen;
        int readStart = (int)rec.qstart;
        int readEnd   = (int)rec.qend;
        char strand   = rec.strand;
        char *path    = rec.path;
        int pathLen   = (int)rec.plen;
        int pathStart = (int)rec.pstart;
        int pathEnd   = (int)rec.pend;
        int qual      = rec.mapq;
        char *cigar   = rec.cigar;

        // Discard reads that have no quality.
        if (qual == 0) continue;
        // A path and CIGAR are required to project the alignment.
        if (!path || !cigar) continue;

        // Ignore reads already mapped elsewhere (first mapping wins).
        khint_t k = mapstr_get(h2, readName);
        if (k < kh_end(h2)) continue;

        // Collect the visited nodes (needed to handle reverse-complement reads).
        uint64_t *nodes;
        int node_size = parse_nodes(readName, path, &nodes, g);
        if (!node_size) continue;

        // If read mapped as rc, reverse the node order and mirror coordinates.
        if (path[0] == '<') {
            reverse(nodes, node_size);
            {
                int new_start = 0, new_end = 0;
                fix_indices(pathStart, pathEnd, pathLen, &new_start, &new_end);
                pathStart = new_start;
                pathEnd = new_end;
            }
            {
                int new_start = 0, new_end = 0;
                fix_indices(readStart, readEnd, readLen, &new_start, &new_end);
                readStart = new_start;
                readEnd = new_end;
            }
        }

        // Initialize the cigar strings.
        char ops[MAX_CIGAR]; int op_index = 0;

        // Soft clip for the prefix part of the read.
        for (int i = 0; i < readStart; i++) ops[op_index++] = __CIGAR_SOFT_CLIP;

        char cigar_ops[MAX_CIGAR];
        int n_ops = parse_cigar(cigar, cigar_ops, MAX_CIGAR, path[0] == '<');
        if (n_ops < 0) {
            fprintf(stderr, "[ERROR] Unable to parse CIGAR %s in read %s\n", cigar, readName);
            free(nodes);
            continue;
        }

        // Initialize the segments (a segment is "reference" iff ref_name != NULL).
        int node_index = 0;
        const gfa_seg_t *seg = gfa_seg_at(g, gfa_idx(g, nodes[node_index]));
        const gfa_seg_t *temp_prev_seg = (seg->ref_name ? seg : NULL);
        int p_length = (int)seg->len - pathStart;

        int reference_start = (seg->ref_name ? seg->start + pathStart + 1 : -1);
        int cigar_op_index = 0, is_reference_start_set = 0;

        while (cigar_op_index < n_ops) {
            char op = cigar_ops[cigar_op_index++];

            if (!is_reference_start_set && seg->ref_name && CIGAR_REF(op)) {
                if (op == __CIGAR_SEQUENCE_MATCH) is_reference_start_set = 1;
                else reference_start++;
            }

            if (seg->ref_name) ops[op_index++] = op;
            else if (CIGAR_QUE(op)) ops[op_index++] = __CIGAR_INSERTION;

            if (CIGAR_REF(op)) {
                p_length--;

                if (p_length) continue;

                node_index++;
                if (node_index == node_size) break;
                seg = gfa_seg_at(g, gfa_idx(g, nodes[node_index]));
                p_length = (int)seg->len;
                if (seg->ref_name) {
                    if (!is_reference_start_set) reference_start = seg->start;

                    if (temp_prev_seg != NULL) {
                        for (int index = temp_prev_seg->end; index < seg->start; index++)
                            ops[op_index++] = __CIGAR_DELETION;
                    }

                    temp_prev_seg = seg;
                }
            }
        }

        // Soft clip for the end part of the read.
        for (int i = readEnd; i < readLen; i++) ops[op_index++] = __CIGAR_SOFT_CLIP;

        #ifdef DEBUG
        // Validate cigar
        int cigar_query_count = 0;
        for (int i = 0; i < op_index; i++) {
            if (CIGAR_QUE(ops[i])) cigar_query_count++;
        }
        if (cigar_query_count != readLen)
            fprintf(stderr, "[ERROR] CIGAR (query) mismatch in length. %d %d\n", cigar_query_count, readLen);
        #endif

        char op = ops[0];
        int op_size = 1;
        for (int i = 1; i < op_index; i++) {
            if (ops[i] == ops[i-1]) {
                op_size++;
            } else {
                if (op == __CIGAR_INSERTION && __SV_MIN_LENGTH <= op_size) {
                    // Insertion detected
                } else if (op == __CIGAR_DELETION && __SV_MIN_LENGTH <= op_size) {
                    // Deletion detected
                } else if (op == __CIGAR_SOFT_CLIP) {
                    // Soft clip at the begining detected
                }
                op = ops[i];
                op_size = 1;
            }
        }
        if (ops[op_index-1] == __CIGAR_SOFT_CLIP && seg->ref_name) {
            // Soft clip at the end detected
        }

        available(gaf_breakpoint, temp_breakpoints, temp_breakpoint_size, temp_breakpoint_capacity);
        temp_breakpoints[temp_breakpoint_size].readStart = readStart;
        temp_breakpoints[temp_breakpoint_size].readEnd = readEnd;
        temp_breakpoints[temp_breakpoint_size].rc = strand == '+' ? 1 : -1;
        temp_breakpoints[temp_breakpoint_size].offset = 0; // TODO: fix offset for correct consensus
        temp_breakpoints[temp_breakpoint_size].type = SV_INS; // TODO: determine SV type

        // Record the read name -> breakpoint index. The key is copied because
        // the reader reuses/frees rec.qname on the next read.
        int absent;
        char *key = strdup(readName);
        validate(key)
        k = mapstr_put(h2, key, &absent);
        if (absent) kh_val(h2, k) = temp_breakpoint_size;
        else free(key);

        free(nodes);
        temp_breakpoint_size++;
    }

    gaf_rec_clear(&rec);
    gaf_close(reader);

    *gaf_breakpoints = temp_breakpoints;
    return temp_breakpoint_size;
}

int parse_fq(const char* file_path, gaf_breakpoint* gaf_breakpoints, mapstr_t *h2) {

    gzFile in = gzopen(file_path, "r");
    if (in == NULL) {
        fprintf(stderr, "Error opening file %s", file_path);
        return 1;
    }

    kseq_t *seq = kseq_init(in);

    while (kseq_read(seq) >= 0) {
        khint64_t k = mapstr_get(h2, seq->name.s);
        if (k == kh_end(h2))
            continue;

        // gaf_breakpoint* aln = gaf_breakpoints + kh_val(h2, k);
        (void)(gaf_breakpoints + kh_val(h2, k));

        // TODO: get substring from the read and prepare for MSA
    }

    kseq_destroy(seq);

    return 0;
}

int discover(int argc, char *argv[]) {

    disc_args params;
    init_disc(argc, argv, &params);

    // Read the graph with links (to validate overlaps) and paths (reference layout).
    gfa_t *g = gfa_read(params.gfa_file, GFA_LINKS | GFA_PATHS);
    if (!g) {
        fprintf(stderr, "[ERROR] GFA file parsing failed.\n");
        return 0;
    }

    // svtrek requires zero-overlap links.
    for (int32_t i = 0; i < gfa_n_link(g); i++) {
        if (gfa_link_at(g, i)->overlap != 0) {
            fprintf(stderr, "[ERROR] Overlaps are not zero, cannot make conversion.\n");
            gfa_destroy(g);
            return 0;
        }
    }

    mapstr_t *h3 = mapstr_init();

    gaf_breakpoint *gaf_breakpoints = NULL;
    // parse_gaf returns the breakpoint count (>= 0), or -1 if the file could
    // not be opened. (The original test here was inverted, aborting on success.)
    int gaf_breakpoint_size = parse_gaf(params.gaf_file, g, &gaf_breakpoints, h3);
    if (gaf_breakpoint_size < 0) {
        fprintf(stderr, "[ERROR] GAF file parsing failed.\n");
        mapstr_destroy(h3);
        gfa_destroy(g);
        return 0;
    }

    parse_fq(params.fq_file, gaf_breakpoints, h3);

    free(gaf_breakpoints);

    // Free the strdup'd read-name keys, then the maps and graph.
    khint_t ki;
    kh_foreach(h3, ki) { free((char *)kh_key(h3, ki)); }
    mapstr_destroy(h3);
    gfa_destroy(g);

    return 1;
}