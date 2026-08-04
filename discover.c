#include "discover.h"

// KHASHL_MAP_INIT(SCOPE, HType, prefix, khkey_t, kh_val_t, __hash_fn, __hash_eq)
KHASHL_MAP_INIT(KH_LOCAL, map32_t, map32, uint64_t, uint32_t, kh_hash_uint64, kh_eq_generic)
KHASHL_MAP_INIT(KH_LOCAL, mapstr_t, mapstr, const char*, uint32_t, kh_hash_str, kh_eq_str)

KSEQ_INIT(gzFile, gzread)

int parse_nodes(const char *readName, char *readPath, uint64_t **nodes, segment *segments, map32_t *h1) {

    uint64_t *temp_nodes = (uint64_t *)malloc(sizeof(uint64_t) * MAX_CIGAR);
    int node_size = 0, path_index, fwd_strand_count = 0, rev_strand_count = 0;
    const char *path = readPath; uint64_t id = 0; char strand = '>';
    
    while ((path_index = next_node(path, &id, &strand))) {
        path += path_index;
        khint_t k = map32_get(h1, id);
        if (k == kh_end(h1)) {
            fprintf(stderr, "[ERROR] Segment %lu in read %s not found.\n", id, readName); 
            free(temp_nodes);
            return 0;
        }
        temp_nodes[node_size++] = id;

        if (strand == '>') fwd_strand_count++;
        else rev_strand_count++;
        // Validation 1: Does node exists
        segment *temp = segments + kh_val(h1, k);
        if (temp->rank > 1) {
            fprintf(stderr, "[ERROR] Read %s contains invalid ranks.\n", readName);
            free(temp_nodes);
            return 0;
        }
        // Validation 2: Does read aligned consistently to one strand
        if (fwd_strand_count && rev_strand_count) {
            fprintf(stderr, "[ERROR] Read %s aligned in mixed strands.\n", readName);
            free(temp_nodes);
            return 0;
        }
    }

    *nodes = temp_nodes;
    return node_size;
}

int parse_gaf(const char* file_path, segment *segments, map32_t *h1, gaf_breakpoint **gaf_breakpoints, mapstr_t *h2) {

    size_t line_cap = MAX_LINE;
    char *line = NULL;
    FILE *file = io_open(file_path, &line, line_cap, "r");
    int line_len = 0;
    
    int temp_breakpoint_size = 0, temp_breakpoint_capacity = 1000000;
    gaf_breakpoint *temp_breakpoints = (gaf_breakpoint *)malloc(temp_breakpoint_capacity * sizeof(gaf_breakpoint));
    validate(temp_breakpoints)

    while ((line_len = io_read(file, &line, &line_cap))) {

        alignment aln;
        memset(&aln, 0, sizeof(alignment));

        char *saveptr, *token;
        int field = 0;

        token = strtok_r(line, "\t", &saveptr);
        while (token != NULL) {
            switch (field) {
                case 0: aln.readName = token; break;
                case 1: aln.readLen = atoi(token); break;
                case 2: aln.readStart = atoi(token); break;
                case 3: aln.readEnd = atoi(token); break;
                case 4: aln.strand = token[0]; break;
                case 5: aln.path = strdup(token); break;
                case 6: aln.pathLen = atoi(token); break;
                case 7: aln.pathStart = atoi(token); break;
                case 8: aln.pathEnd = atoi(token); break;
                case 9: aln.matches = atoi(token); break;
                case 10: aln.blockLen = atoi(token); break;
                case 11: aln.qual = atoi(token); break;
                default:
                    if (strncmp(token, "cg:Z:", 5) == 0) aln.cigar = strdup(token + 5);
                    break;
            }
            field++;
            token = strtok_r(NULL, "\t", &saveptr);
        }

        // Discard reads that has no quality
        if (aln.qual == 0) {
            free(aln.path);
            if (aln.cigar) free(aln.cigar);
            continue;
        }

        // Check if read is already mapped to somewhere else. If yes, ignore the alignment
        // TODO: check if this mapping has more qual score or smt else so it can be replaced
        khint_t k = mapstr_get(h2, aln.readName);
        if (k < kh_end(h2)) {
            free(aln.path);
            if (aln.cigar) free(aln.cigar);
            continue;
        }

        // Get nodes and store in array. It is needed in case if it is mapped as reverse complement
        uint64_t *nodes;
        int node_size = parse_nodes(aln.readName, aln.path, &nodes, segments, h1);
        // Continue if the wile loop terminated because of missing segment
        if (!node_size) {
            free(aln.path);
            if (aln.cigar) free(aln.cigar);
            continue;
        }

        // If read mapped as rc, reverse the order of the nodes
        if (aln.path[0] == '<') {
            reverse(nodes, node_size);
            {
                int new_start = 0, new_end = 0;
                fix_indices(aln.pathStart, aln.pathEnd, aln.pathLen, &new_start, &new_end);
                aln.pathStart = new_start;
                aln.pathEnd = new_end;
            }
            {
                int new_start = 0, new_end = 0;
                fix_indices(aln.readStart, aln.readEnd, aln.readLen, &new_start, &new_end);
                aln.readStart = new_start;
                aln.readEnd = new_end;
            }
        }

        // Initialize the cigar strings
        char ops[MAX_CIGAR]; int op_index = 0;

        // Soft clip for the prefix part of the read
        for (int i = 0; i < aln.readStart; i++) ops[op_index++] = __CIGAR_SOFT_CLIP;

        char cigar_ops[MAX_CIGAR];
        int n_ops = parse_cigar(aln.cigar, cigar_ops, MAX_CIGAR, aln.path[0] == '<');
        if (n_ops < 0) {
            fprintf(stderr, "[ERROR] Unable to parse CIGAR %s in read %s\n", aln.cigar, aln.readName);
            free(aln.path); free(nodes);
            if (aln.cigar) free(aln.cigar);
            continue;
        }

        // Initialize the segments
        int node_index = 0;
        k = map32_get(h1, nodes[node_index]);
        segment *seg = segments + kh_val(h1, k);
        segment *temp_prev_seg = (seg->rank == 0 ? seg : NULL);
        int p_length = strlen(seg->seq)-aln.pathStart;

        int reference_start = (seg->rank == 0 ? seg->start+aln.pathStart+1 : -1);
        int cigar_op_index = 0, is_reference_start_set = 0;

        while (cigar_op_index < n_ops) {
            char op = cigar_ops[cigar_op_index++];

            if (!is_reference_start_set && seg->rank == 0 && CIGAR_REF(op)) {
                if (op == __CIGAR_SEQUENCE_MATCH) is_reference_start_set = 1;
                else reference_start++;
            }

            if (seg->rank == 0) ops[op_index++] = op;
            else if (CIGAR_QUE(op)) ops[op_index++] = __CIGAR_INSERTION;

            if (CIGAR_REF(op)) {
                p_length--;

                if (p_length) continue;

                node_index++;
                if (node_index == node_size) break;
                khint_t k = map32_get(h1, nodes[node_index]);
                seg = segments + kh_val(h1, k);
                p_length = strlen(seg->seq);
                if (seg->rank == 0) {
                    if(!is_reference_start_set) reference_start = seg->start;

                    if (temp_prev_seg != NULL) {
                        for (int index = temp_prev_seg->end; index < seg->start; index++)
                            ops[op_index++] = __CIGAR_DELETION;
                    }

                    temp_prev_seg = seg;
                }
            }
        }

        // Soft clip for the end part of the read
        for (int i = aln.readEnd; i < aln.readLen; i++) ops[op_index++] = __CIGAR_SOFT_CLIP;

        #ifdef DEBUG
        // Validate cigar
        int cigar_query_count = 0;
        for (int i=0; i<op_index; i++) {
            if (CIGAR_QUE(ops[i])) cigar_query_count++;
        }
        if (cigar_query_count != aln.readLen) 
            fprintf(stderr, "[ERROR] CIGAR (query) mismatch in length. %d %d\n", cigar_query_count, aln.readLen);
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
        if (ops[op_index-1] == __CIGAR_SOFT_CLIP && seg->rank == 0) {
            // Soft clip at the end detected
        }

        // Cleanup
        free(aln.path); free(nodes);
        if (aln.cigar) free(aln.cigar);

        available(gaf_breakpoint, temp_breakpoints, temp_breakpoint_size, temp_breakpoint_capacity);
        temp_breakpoints[temp_breakpoint_size].readStart = aln.readStart;
        temp_breakpoints[temp_breakpoint_size].readEnd = aln.readEnd;
        temp_breakpoints[temp_breakpoint_size].rc = aln.strand == '>' ? 1 : -1;
        temp_breakpoints[temp_breakpoint_size].offset = 0; // TODO: fix offset for correct consensus
        temp_breakpoints[temp_breakpoint_size].type = SV_INS; // TODO: determine SV type

        int absent;
        k = mapstr_put(h2, aln.readName, &absent);
        kh_val(h2, k) = temp_breakpoint_size;

        temp_breakpoint_size++;
    }    
    fclose(file);
    free(line);

    *gaf_breakpoints = temp_breakpoints;
    return temp_breakpoint_size;
}

int parse_gfa(const char* file_path, segment **segments, int *segment_size, map32_t *h1) {

    // NOTE: Rank for each segment is set to be 1 as default.
    // NOTE: It is assumed that ranks can be only 0 and 1. No rank > 1 is allowed
    // NOTE: No overlap > 0 is allowed

    size_t line_cap = MAX_LINE;
    char *line = NULL;
    FILE *file = io_open(file_path, &line, line_cap, "r");
    int line_len = 0;
    
    // Allocation for Segment
    int temp_segment_size = 0, segment_capacity = 1000000;
    segment *temp_segments = (segment *)malloc(segment_capacity * sizeof(segment));
    validate(temp_segments)

    while ((line_len = io_read(file, &line, &line_cap))) {
        if (line[0] == 'S') {
            
            available(segment, temp_segments, temp_segment_size, segment_capacity)
            
            char *saveptr;
            char *token = strtok_r(line, "\t", &saveptr); // Skip 'S'
            // ID
            token = strtok_r(NULL, "\t", &saveptr); // ID
            temp_segments[temp_segment_size].id = strtoull(token, NULL, 10);
            
            // Sequence
            token = strtok_r(NULL, "\t", &saveptr); // Sequence string
            size_t seq_len = strlen(token);
            temp_segments[temp_segment_size].seq = malloc(seq_len + 1);
            validate(temp_segments[temp_segment_size].seq)
            strcpy(temp_segments[temp_segment_size].seq, token);
            
            temp_segments[temp_segment_size].rank = 1;
            temp_segments[temp_segment_size].start = -1; 
            temp_segments[temp_segment_size].end = seq_len;
            temp_segments[temp_segment_size].next = NULL;
            
            khint_t k; int absent;
            k = map32_put(h1, temp_segments[temp_segment_size].id, &absent);
            kh_val(h1, k) = temp_segment_size; // set value as index in segments
            temp_segment_size++;
        } else if (line[0] == 'L') {
            // Do nothing at this point
        } else if (line[0] == 'P') {
            char *token = strtok(line, "\t");
            token = strtok(NULL, "\t"); // process path ID
            token = strtok(NULL, "\t"); // process segment IDs

            char *segment_token = strtok(token, ",");
            int ref_pos = 0;

            while (segment_token) {
                size_t len = strlen(segment_token);
                if (segment_token[len - 1] == '+') segment_token[len - 1] = '\0'; // remove the trailing '+'

                uint64_t seg_id = strtoull(segment_token, NULL, 10);
                khint_t k = map32_get(h1, seg_id);
                segment *current_segment = temp_segments + kh_val(h1, k);
                
                current_segment->rank = 0;
                current_segment->start = ref_pos;
                ref_pos += strlen(current_segment->seq);
                current_segment->end = ref_pos;
                segment_token = strtok(NULL, ",");  // move to next segment
            }
        }
    }

    rewind(file);

    // Set start pos of alternative paths' nodes. It is needed for MSE for large INS/ALT/DUP
    // It is assumed that the links are ordered w.r.t. original and relative positions
    while ((line_len = io_read(file, &line, &line_cap))) {
        if (line[0] == 'L') {
            uint64_t id1, id2;
            char strand1, strand2;
            size_t overlap = 0;
            sscanf(line, "L\t%lu\t%c\t%lu\t%c\t%luM", &id1, &strand1, &id2, &strand2, &overlap);

            // Here, we validate overlaps, making sure that they are none (0)
            if (overlap) {
                fprintf(stderr, "[ERROR] Overlaps are not zero, cannot make conversion.\n");
                exit(1);
            }

            khint_t k1 = map32_get(h1, id1);
            if (kh_val(h1, k1) == kh_end(h1)) {
                fprintf(stderr, "[ERROR] Segment %lu does not exists.\n", id1);
                exit(1);
            }
            segment *segment1 = temp_segments + kh_val(h1, k1);

            khint_t k2 = map32_get(h1, id2);
            if (kh_val(h1, k2) == kh_end(h1)) {
                fprintf(stderr, "[ERROR] Segment %lu does not exists.\n", id2);
                exit(1);
            }
            segment *segment2 = temp_segments + kh_val(h1, k2);

            if (segment1->rank && segment2->rank) {
                segment1->next = segment2;
            } else if (segment1->rank == 0 && segment2->rank) {
                segment2->start = 0;
            }
        }
    }

    io_close(file, &line);

    for (int i = 0; i < temp_segment_size; i++) {
        if (temp_segments[i].rank == 1 && temp_segments[i].start == 0) {
            int path_length = 0;

            segment *current_segment = temp_segments + i;

            while (current_segment != NULL) {
                current_segment->start = path_length;
                path_length += current_segment->end;
                current_segment->end = path_length;
                current_segment = current_segment->next;
            }
        }
    }

    *segments = temp_segments;
    *segment_size = temp_segment_size;

    return 0;
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

    map32_t *h1 = map32_init();

    segment *segments; int segment_size;
    if (parse_gfa(params.gfa_file, &segments, &segment_size, h1)) {
        fprintf(stderr, "[ERROR] GFA file parsing resulted in null.\n");
        return 0;
    }

    mapstr_t *h3 = mapstr_init();

    gaf_breakpoint *gaf_breakpoints; int gaf_breakpoint_size;
    if ((gaf_breakpoint_size = parse_gaf(params.gaf_file, segments, h1, &gaf_breakpoints, h3))) {
        fprintf(stderr, "[ERROR] GAF file parsing failed.\n");
        return 0;
    }

    parse_fq(params.fq_file, gaf_breakpoints, h3);

    // Cleanup segments
    for (int i = 0; i < segment_size; i++) {
        if (segments[i].seq) free(segments[i].seq);
    }
    free(segments);
    free(gaf_breakpoints);

    map32_destroy(h1);
    mapstr_destroy(h3);

    return 1;
}