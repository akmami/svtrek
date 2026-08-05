#include "refinement.h"
#include "msa.h"
#include <htslib/sam.h>


#define abs_diff(a, b) (a < b ? (b - a) : (a - b))


/* ------------------------------------------------------------------ */
/*  Consensus of a set of candidate breakpoint locations              */
/* ------------------------------------------------------------------ */

/*
 * Legacy, currently-unused helper kept for reference. Finds the value that
 * anchors the largest cluster (values within `consensus_interval`).
 */
svtrek_index consensus(svtrek_index *arr, size_t size, size_t consensus_min_count, svtrek_index consensus_interval) {
    svtrek_index consensus_val = 0;
    size_t max_count = consensus_min_count - 1;

    quicksort(arr, 0, size-1);

    for (size_t i = 0; i < size; i++) {
        size_t count = 1;
        // count how many values fall within consensus_interval of arr[i]
        for (size_t j = i + 1; j < size && arr[j] <= arr[i] + consensus_interval; j++) {
            count++;
        }
        if (count > max_count) {
            max_count = count;
            consensus_val = arr[i];
        }
    }
    return consensus_val;
}

/*
 * Given a set of candidate breakpoint `locations` and an imprecise reference
 * position `pos`, return a refined position.
 *
 * A candidate is the (rounded) mean of a cluster of locations that all fall
 * within `consensus_interval` of each other. Only clusters whose anchor is
 * within `consensus_interval_range` of `pos` and that reach at least
 * `consensus_min_count` supporting reads are considered. Among the eligible
 * clusters we prefer the one CLOSEST to `pos`, breaking ties by higher
 * support -- this keeps refinement anchored to the annotated location
 * instead of snapping onto a stronger but unrelated neighbouring variant.
 * Returns -1 when nothing qualifies.
 *
 * NOTE: this replaces an earlier version that relied on two buggy
 * lower_bound/upper_bound helpers (upper_bound returned index 0 for almost
 * every input, so the right-hand scan effectively never ran).
 */
svtrek_index consensus_pos(svtrek_index *locations, size_t size, svtrek_index pos, size_t consensus_min_count, svtrek_index consensus_interval, svtrek_index consensus_interval_range) {

    if (size < consensus_min_count) {
        return -1;
    }

    quicksort(locations, 0, size - 1);

    svtrek_index best_candidate = 0;
    svtrek_index best_distance = 0x7fffffff;
    size_t best_count = 0;

    for (size_t i = 0; i < size; i++) {
        /* Grow a cluster of locations within consensus_interval of one
         * another, anchored at locations[i] (array is sorted). */
        size_t count = 1;
        int64_t total = locations[i];
        for (size_t j = i + 1; j < size && locations[j] <= locations[i] + consensus_interval; j++) {
            count++;
            total += locations[j];
        }

        if (count < consensus_min_count)
            continue;

        svtrek_index candidate = (svtrek_index)((total + count / 2) / count);
        svtrek_index distance = abs_diff(pos, candidate);

        /* Reject clusters whose refined position is farther than the allowed
         * relocation window from the imprecise (annotated) position. */
        if (distance >= consensus_interval_range)
            continue;

        /* Stick to the exact location when possible: prefer the cluster
         * CLOSEST to the annotated position, using support only to break
         * ties. This stops refinement from snapping onto a stronger but
         * unrelated neighbouring variant that merely falls inside the window. */
        if (distance < best_distance || (distance == best_distance && count > best_count)) {
            best_distance = distance;
            best_count = count;
            best_candidate = candidate;
        }
    }

    return best_candidate;
}

/* ------------------------------------------------------------------ */
/*  Breakpoint-position refinement (CIGAR based)                      */
/* ------------------------------------------------------------------ */

svtrek_index refine_start(sv_type_t sv_type, int chrom, interval inter, uint32_t imprecise_pos, t_arg *params) {

    size_t capacity = 100;
    size_t size = 0;
    svtrek_index *start_locations = (svtrek_index *)malloc(sizeof(svtrek_index)*capacity);
    if(start_locations == NULL) {
        fprintf(stderr, "Couldn't allocate array for start positions.\n");
        return -1;
    }

    bam1_t *aln = bam_init1();
    int tid = chrom;
    hts_itr_t *iter = sam_itr_queryi(params->hargs.bam_file_index, tid, inter.start-1, inter.end-1);

    if (iter) {
        while (sam_itr_next(params->hargs.fp_in, iter, aln) > 0) {
            uint32_t reference_pos = aln->core.pos;
            uint32_t *cigar = bam_get_cigar(aln);
            int check_soft_clip = bam_cigar_op(cigar[aln->core.n_cigar-1]) == __CIGAR_SOFT_CLIP;

            // n_cigar is unsigned, so use uint32_t in the loop
            for (uint32_t i = 0; i < aln->core.n_cigar; i++) {
                if (sv_type == SV_DEL && bam_cigar_op(cigar[i]) == __CIGAR_DELETION && __SV_MIN_LENGTH < bam_cigar_oplen(cigar[i])) {
                    if (capacity == size) {
                        capacity = capacity * 1.5;
                        svtrek_index *temp = (svtrek_index *)realloc(start_locations, sizeof(int)*capacity);
                        if (temp == NULL) {
                            fprintf(stderr, "[ERROR] Couldn't reallocate start locations array.\n");
                            free(start_locations);
                            bam_destroy1(aln);
                            sam_itr_destroy(iter);
                            return -1;
                        }
                        start_locations = temp;
                    }
                    // reference_pos here is the 0-based first deleted base,
                    // which numerically equals the 1-based anchor (VCF POS).
                    start_locations[size++] = reference_pos;
                }

                if (bam_cigar_op(cigar[i]) != __CIGAR_INSERTION && bam_cigar_op(cigar[i]) != __CIGAR_SOFT_CLIP) {
                    reference_pos += bam_cigar_oplen(cigar[i]);
                }

                if (reference_pos > inter.end) {
                    check_soft_clip = 0;
                    break;
                }
            }

            if (check_soft_clip && inter.start <= reference_pos && reference_pos <= inter.end) {
                if (capacity == size) {
                    capacity = capacity * 1.5;
                    svtrek_index *temp = (svtrek_index *)realloc(start_locations, sizeof(svtrek_index)*capacity);
                    if (temp == NULL) {
                        fprintf(stderr, "[ERROR] Couldn't reallocate start locations array.\n");
                        free(start_locations);
                        bam_destroy1(aln);
                        sam_itr_destroy(iter);
                        return -1;
                    }
                    start_locations = temp;
                }
                // A read soft-clipped at its 3' end stops aligning at the
                // deletion start; reference_pos is that (0-based) coordinate.
                start_locations[size++] = reference_pos;
            }
        }
        sam_itr_destroy(iter);
    }
    bam_destroy1(aln);

    // pick the best refined start
    svtrek_index result = consensus_pos(start_locations, size, imprecise_pos, params->consensus_min_count, params->consensus_interval, params->consensus_interval_range);
    free(start_locations);
    return result;
}

svtrek_index refine_end(sv_type_t sv_type, int chrom, interval inter, uint32_t imprecise_pos, t_arg *params) {

    size_t capacity = 100;
    size_t size = 0;
    svtrek_index *end_locations = (svtrek_index *)malloc(sizeof(svtrek_index)*capacity);
    if(end_locations == NULL) {
        fprintf(stderr, "Couldn't allocate array for end positions.\n");
        return -1;
    }

    bam1_t *aln = bam_init1();
    int tid = chrom;
    hts_itr_t *iter = sam_itr_queryi(params->hargs.bam_file_index, tid, inter.start-1, inter.end-1);

    if (iter) {
        while (sam_itr_next(params->hargs.fp_in, iter, aln) > 0) {
            uint32_t reference_pos = aln->core.pos;
            uint32_t *cigar = bam_get_cigar(aln);

            for (uint32_t i = 0; i < aln->core.n_cigar; i++) {
                if (sv_type == SV_DEL && bam_cigar_op(cigar[i]) == __CIGAR_DELETION && __SV_MIN_LENGTH < bam_cigar_oplen(cigar[i])) {
                    if (capacity == size) {
                        capacity = capacity * 1.5;
                        svtrek_index *temp = (svtrek_index *)realloc(end_locations, sizeof(svtrek_index)*capacity);
                        if (temp == NULL) {
                            fprintf(stderr, "[ERROR] Couldn't reallocate end locations array.\n");
                            free(end_locations);
                            bam_destroy1(aln);
                            sam_itr_destroy(iter);
                            return -1;
                        }
                        end_locations = temp;
                    }
                    // 0-based first deleted base + deletion length == the
                    // 1-based last deleted base (VCF END). The previous "+1"
                    // over-shot the end by one base.
                    end_locations[size++] = reference_pos + bam_cigar_oplen(cigar[i]);
                }

                if (bam_cigar_op(cigar[i]) != __CIGAR_INSERTION && bam_cigar_op(cigar[i]) != __CIGAR_SOFT_CLIP) {
                    reference_pos += bam_cigar_oplen(cigar[i]);
                }

                if (reference_pos > inter.end) {
                    break;
                }
            }

            if (bam_cigar_op(cigar[0]) == __CIGAR_SOFT_CLIP && inter.start <= aln->core.pos && (uint32_t)aln->core.pos <= inter.end) {
                if (capacity == size) {
                    capacity = capacity * 1.5;
                    svtrek_index *temp = (svtrek_index *)realloc(end_locations, sizeof(svtrek_index)*capacity);
                    if (temp == NULL) {
                        fprintf(stderr, "[ERROR] Couldn't reallocate end locations array.\n");
                        free(end_locations);
                        bam_destroy1(aln);
                        sam_itr_destroy(iter);
                        return -1;
                    }
                    end_locations = temp;
                }
                // A read soft-clipped at its 5' end resumes aligning right
                // after the deletion; aln->core.pos is that coordinate
                // (0-based), numerically the 1-based last deleted base.
                // The previous code stored the alignment END position here,
                // which is unrelated to the breakpoint.
                end_locations[size++] = aln->core.pos;
            }
        }
        sam_itr_destroy(iter);
    }
    bam_destroy1(aln);

    // pick the best refined end
    svtrek_index result = consensus_pos(end_locations, size, imprecise_pos, params->consensus_min_count, params->consensus_interval, params->consensus_interval_range);
    free(end_locations);
    return result;
}

svtrek_index refine_point(sv_type_t sv_type, int chrom, interval inter, svtrek_index imprecise_pos, t_arg *params) {

    size_t capacity = 100;
    size_t size = 0;
    svtrek_index *locations = (svtrek_index *)malloc(sizeof(svtrek_index)*capacity);
    if(locations == NULL) {
        fprintf(stderr, "Couldn't allocate array for positions.\n");
        return -1;
    }

    bam1_t *aln = bam_init1();

    int tid = chrom;
    hts_itr_t *iter = sam_itr_queryi(params->hargs.bam_file_index, tid, inter.start-1, inter.end-1);
    if (iter) {
        while (sam_itr_next(params->hargs.fp_in, iter, aln) > 0) {
            uint32_t reference_pos = aln->core.pos;
            uint32_t *cigar   = bam_get_cigar(aln);

            for (uint32_t i = 0; i < aln->core.n_cigar; i++) {
                // NOTE: robust inversion detection needs supplementary-alignment
                // (SA tag) strand analysis and is not implemented yet. As a
                // best-effort signal we record large soft-clips, which mark the
                // breakpoints of a split read.
                if (sv_type == SV_INV && bam_cigar_op(cigar[i]) == __CIGAR_SOFT_CLIP && __SV_MIN_LENGTH < bam_cigar_oplen(cigar[i])) {
                    if (capacity == size) {
                        capacity = capacity * 1.5;
                        svtrek_index *temp = (svtrek_index *)realloc(locations, sizeof(svtrek_index)*capacity);
                        if (temp == NULL) {
                            fprintf(stderr, "[ERROR] Couldn't reallocate locations array.\n");
                            free(locations);
                            bam_destroy1(aln);
                            sam_itr_destroy(iter);
                            return -1;
                        }
                        locations = temp;
                    }
                    locations[size++] = reference_pos;
                }

                if (bam_cigar_op(cigar[i]) != __CIGAR_INSERTION && bam_cigar_op(cigar[i]) != __CIGAR_SOFT_CLIP) {
                    reference_pos += bam_cigar_oplen(cigar[i]);
                }

                if (reference_pos > inter.end) {
                    break;
                }
            }
        }
        sam_itr_destroy(iter);
    }
    bam_destroy1(aln);
    svtrek_index result = consensus_pos(locations, size, imprecise_pos, params->consensus_min_count, params->consensus_interval, params->consensus_interval_range);
    free(locations);
    return result;
}

svtrek_index refine_ins(int chrom, interval inter, svtrek_index imprecise_pos, t_arg *params) {

    size_t capacity = 100;
    size_t size = 0;
    svtrek_index *locations = (svtrek_index *)malloc(sizeof(svtrek_index)*capacity);
    if(locations == NULL) {
        fprintf(stderr, "Couldn't allocate array for positions.\n");
        return -1;
    }

    bam1_t *aln = bam_init1();

    int tid = chrom;
    hts_itr_t *iter = sam_itr_queryi(params->hargs.bam_file_index, tid, inter.start-1, inter.end-1);
    if (iter) {
        while (sam_itr_next(params->hargs.fp_in, iter, aln) > 0) {
            uint32_t reference_pos = aln->core.pos;
            uint32_t *cigar   = bam_get_cigar(aln);

            for (uint32_t i = 0; i < aln->core.n_cigar; i++) {
                if (bam_cigar_op(cigar[i]) == __CIGAR_INSERTION && __SV_MIN_LENGTH <= bam_cigar_oplen(cigar[i])) {
                    if (capacity == size) {
                        capacity = capacity * 1.5;
                        svtrek_index *temp = (svtrek_index *)realloc(locations, sizeof(svtrek_index)*capacity);
                        if (temp == NULL) {
                            fprintf(stderr, "[ERROR] Couldn't reallocate locations array.\n");
                            free(locations);
                            bam_destroy1(aln);
                            sam_itr_destroy(iter);
                            return -1;
                        }
                        locations = temp;
                    }
                    // reference_pos at an insertion is the anchor base
                    // (numerically the 1-based VCF POS).
                    locations[size++] = reference_pos;
                }

                if (bam_cigar_op(cigar[i]) != __CIGAR_INSERTION && bam_cigar_op(cigar[i]) != __CIGAR_SOFT_CLIP) {
                    reference_pos += bam_cigar_oplen(cigar[i]);
                }

                if (reference_pos > inter.end) {
                    break;
                }
            }
        }
        sam_itr_destroy(iter);
    }
    bam_destroy1(aln);
    svtrek_index result = consensus_pos(locations, size, imprecise_pos, params->consensus_min_count, params->consensus_interval, params->consensus_interval_range);
    free(locations);
    return result;
}

/* ------------------------------------------------------------------ */
/*  Variant-sequence extraction for abPOA MSA                         */
/* ------------------------------------------------------------------ */

/*
 * Collect the inserted query bases (I-op of length >= __SV_MIN_LENGTH) from
 * every read whose insertion anchor falls within `tol` of `cons_pos`. These
 * are the candidate ALT (inserted) sequences that feed the MSA.
 */
static size_t collect_insertion_seqs(int chrom, svtrek_index cons_pos, svtrek_index tol, t_arg *params,
                                     char ***seqs_out, int **lens_out, size_t max_seqs) {
    char **seqs = (char **)malloc(sizeof(char *) * max_seqs);
    int *lens = (int *)malloc(sizeof(int) * max_seqs);
    size_t n = 0;

    if (seqs == NULL || lens == NULL) {
        free(seqs); free(lens);
        *seqs_out = NULL; *lens_out = NULL;
        return 0;
    }

    svtrek_index lo = cons_pos - tol; if (lo < 1) lo = 1;
    svtrek_index hi = cons_pos + tol;

    int tid = chrom;
    bam1_t *aln = bam_init1();
    hts_itr_t *iter = sam_itr_queryi(params->hargs.bam_file_index, tid, lo - 1, hi);
    if (iter) {
        while (n < max_seqs && sam_itr_next(params->hargs.fp_in, iter, aln) > 0) {
            uint32_t ref = aln->core.pos;
            uint32_t qpos = 0;
            uint32_t *cigar = bam_get_cigar(aln);
            uint8_t  *qseq  = bam_get_seq(aln);

            for (uint32_t i = 0; i < aln->core.n_cigar && n < max_seqs; i++) {
                int op  = bam_cigar_op(cigar[i]);
                size_t len = bam_cigar_oplen(cigar[i]);

                if (op == __CIGAR_INSERTION && len >= __SV_MIN_LENGTH && abs_diff(ref, cons_pos) <= tol) {
                    char *s = (char *)malloc(len + 1);
                    if (s != NULL) {
                        for (size_t j = 0; j < len; j++)
                            s[j] = seq_nt16_str[bam_seqi(qseq, qpos + j)];
                        s[len] = '\0';
                        seqs[n] = s;
                        lens[n] = (int)len;
                        n++;
                    }
                }

                if (bam_cigar_type(op) & 1) qpos += len;   /* consumes query */
                if (bam_cigar_type(op) & 2) ref  += len;   /* consumes reference */

                if (ref > hi + len) break;
            }
        }
        sam_itr_destroy(iter);
    }
    bam_destroy1(aln);

    *seqs_out = seqs;
    *lens_out = lens;
    return n;
}

/*
 * Collect, for every read that fully spans the reference window
 * [lo-flank, hi+flank), the query subsequence aligned to that window
 * (inserted bases included, deleted reference skipped). For a true deletion
 * these spanning reads jump across the deleted reference, so their consensus
 * length is roughly the flanking length -- shorter than the window by the
 * deletion size.
 */
static size_t collect_spanning_seqs(int chrom, svtrek_index lo, svtrek_index hi, svtrek_index flank, t_arg *params,
                                    char ***seqs_out, int **lens_out, size_t max_seqs) {
    svtrek_index win_lo; 
    if (lo < flank) win_lo = 0;
    else            win_lo = lo - flank;
    svtrek_index win_hi = hi + flank;

    char **seqs = (char **)malloc(sizeof(char *) * max_seqs);
    int *lens = (int *)malloc(sizeof(int) * max_seqs);
    size_t n = 0;

    if (seqs == NULL || lens == NULL) {
        free(seqs); free(lens);
        *seqs_out = NULL; *lens_out = NULL;
        return 0;
    }

    size_t cap_seq = (win_hi - win_lo) + 16;
    int tid = chrom;
    bam1_t *aln = bam_init1();
    hts_itr_t *iter = sam_itr_queryi(params->hargs.bam_file_index, tid, win_lo, win_hi);
    if (iter) {
        while (n < max_seqs && sam_itr_next(params->hargs.fp_in, iter, aln) > 0) {
            svtrek_index aln_start = aln->core.pos;
            svtrek_index aln_end   = bam_endpos(aln);           /* 0-based, exclusive */
            /* require the read to fully span the window */
            if (aln_start > win_lo || aln_end < win_hi)
                continue;

            svtrek_index ref = aln->core.pos;
            uint32_t qpos = 0;
            uint32_t *cigar = bam_get_cigar(aln);
            uint8_t  *qseq  = bam_get_seq(aln);

            char *s = (char *)malloc(cap_seq + 1);
            if (s == NULL) continue;
            size_t sl = 0;

            for (uint32_t i = 0; i < aln->core.n_cigar; i++) {
                int op  = bam_cigar_op(cigar[i]);
                size_t len = bam_cigar_oplen(cigar[i]);
                int cq  = bam_cigar_type(op) & 1;      /* consumes query */
                int cr  = bam_cigar_type(op) & 2;      /* consumes reference */

                if (cq && cr) {                        /* M / = / X */
                    for (size_t j = 0; j < len && sl < cap_seq; j++) {
                        svtrek_index rp = ref + j;
                        if (rp >= win_lo && rp < win_hi)
                            s[sl++] = seq_nt16_str[bam_seqi(qseq, qpos + j)];
                    }
                } else if (op == __CIGAR_INSERTION) {  /* inserted bases anchored at ref */
                    if (ref >= win_lo && ref < win_hi)
                        for (size_t j = 0; j < len && sl < cap_seq; j++)
                            s[sl++] = seq_nt16_str[bam_seqi(qseq, qpos + j)];
                }

                if (cq) qpos += len;
                if (cr) ref  += len;
            }

            if (sl > 0) {
                s[sl] = '\0';
                seqs[n] = s;
                lens[n] = (int)sl;
                n++;
            } else {
                free(s);
            }
        }
        sam_itr_destroy(iter);
    }
    bam_destroy1(aln);

    *seqs_out = seqs;
    *lens_out = lens;
    return n;
}

/* ------------------------------------------------------------------ */
/*  Deletion breakpoint signatures (joint start,end) + clustering      */
/* ------------------------------------------------------------------ */

static int push_index(svtrek_index **arr, size_t *size, size_t *cap, svtrek_index value) {
    if (*size == *cap) {
        size_t new_cap = ((*cap) * 1.5) + 1;
        svtrek_index *tmp = (svtrek_index *)realloc(*arr, sizeof(svtrek_index) * new_cap);
        if (tmp == NULL)
            return -1;
        *arr = tmp;
        *cap = new_cap;
    }
    (*arr)[(*size)++] = value;
    return 0;
}

/*
 * Per-read deletion signatures collected over a search window.
 *
 * A deletion is a single event that each read observes at most partially:
 *   - a read that aligns across it carries a CIGAR D-op, giving BOTH
 *     breakpoints at once  -> a joint (tup_start, tup_end) tuple;
 *   - a read clipped at the junction is a split read, giving ONE side:
 *       * a large 3' (right) soft-clip ends at the deletion start  -> sc_start
 *       * a large 5' (left)  soft-clip resumes at the deletion end  -> sc_end
 *
 * Clustering the tuples and taking the majority keeps start and end coupled;
 * the one-sided soft-clip signals then sharpen each breakpoint (and, when the
 * deletion is longer than a read and no D-op exists, stand in for the tuples).
 */
typedef struct {
    svtrek_index *tup_start, *tup_end;  size_t n_tuple;    /* D-op reads: joint (start,end) */
    svtrek_index *sc_start;             size_t n_sc_start; /* right soft-clip -> start only  */
    svtrek_index *sc_end;               size_t n_sc_end;   /* left  soft-clip -> end   only  */
} del_signatures;

static void free_del_signatures(del_signatures *sg) {
    free(sg->tup_start); free(sg->tup_end);
    free(sg->sc_start);  free(sg->sc_end);
    sg->tup_start = sg->tup_end = sg->sc_start = sg->sc_end = NULL;
    sg->n_tuple = sg->n_sc_start = sg->n_sc_end = 0;
}

/* Append one (start,end) tuple, growing both arrays together so their
 * capacities never drift apart. Returns -1 on OOM. */
static int push_tuple(del_signatures *sg, size_t *cap, svtrek_index start, svtrek_index end) {
    if (sg->n_tuple == *cap) {
        size_t new_cap = (size_t)((*cap) * 1.5) + 1;
        svtrek_index *ts = (svtrek_index *)realloc(sg->tup_start, sizeof(svtrek_index) * new_cap);
        svtrek_index *te = (svtrek_index *)realloc(sg->tup_end,   sizeof(svtrek_index) * new_cap);
        if (ts) sg->tup_start = ts;
        if (te) sg->tup_end   = te;
        if (ts == NULL || te == NULL)
            return -1;
        *cap = new_cap;
    }
    sg->tup_start[sg->n_tuple] = start;
    sg->tup_end[sg->n_tuple]   = end;
    sg->n_tuple++;
    return 0;
}

static del_signatures collect_deletion_signatures(int chrom, interval win, interval annotated, t_arg *params) {

    del_signatures sg;
    size_t cap_tuple = 64, cap_sc_start = 64, cap_sc_end = 64;
    sg.tup_start = (svtrek_index *)malloc(sizeof(svtrek_index) * cap_tuple);
    sg.tup_end   = (svtrek_index *)malloc(sizeof(svtrek_index) * cap_tuple);
    sg.sc_start  = (svtrek_index *)malloc(sizeof(svtrek_index) * cap_sc_start);
    sg.sc_end    = (svtrek_index *)malloc(sizeof(svtrek_index) * cap_sc_end);
    sg.n_tuple = sg.n_sc_start = sg.n_sc_end = 0;

    if (sg.tup_start == NULL || sg.tup_end == NULL || sg.sc_start == NULL || sg.sc_end == NULL) {
        fprintf(stderr, "[ERROR] Couldn't allocate deletion signature arrays.\n");
        free_del_signatures(&sg);
        return sg;
    }

    svtrek_index range = params->consensus_interval_range;

    bam1_t *aln = bam_init1();
    int tid = chrom;
    hts_itr_t *iter = sam_itr_queryi(params->hargs.bam_file_index, tid, (int)win.start - 1, (int)win.end);
    if (iter) {
        while (sam_itr_next(params->hargs.fp_in, iter, aln) > 0) {
            uint32_t *cigar   = bam_get_cigar(aln);
            uint32_t  n_cigar = aln->core.n_cigar;
            if (n_cigar == 0)
                continue;

            // joint tuple: the D-op closest to the annotated event, so one read contributes at most one (start,end) vote
            uint32_t reference_pos = aln->core.pos;
            int have_tuple = 0;
            svtrek_index tuple_start = 0, tuple_end = 0;
            int64_t best_dist = 0x7fffffffffffffffLL;

            for (uint32_t i = 0; i < n_cigar; i++) {
                int op  = bam_cigar_op(cigar[i]);
                size_t len = bam_cigar_oplen(cigar[i]);

                if (op == __CIGAR_DELETION && __SV_MIN_LENGTH < len) {
                    svtrek_index d_start = reference_pos;
                    svtrek_index d_end   = reference_pos + len;
                    int64_t dist = (int64_t)abs_diff(annotated.start, d_start) + (int64_t)abs_diff(annotated.end, d_end);
                    if (dist < best_dist) {
                        best_dist   = dist;
                        tuple_start = d_start;
                        tuple_end   = d_end;
                        have_tuple  = 1;
                    }
                }

                if (op != __CIGAR_INSERTION && op != __CIGAR_SOFT_CLIP)
                    reference_pos += len;

                if (reference_pos > win.end)
                    break;
            }

            if (have_tuple)
                if (push_tuple(&sg, &cap_tuple, tuple_start, tuple_end) < 0)
                    break;

            // right (3') soft-clip; deletion START at the alignment end
            if (bam_cigar_op(cigar[n_cigar - 1]) == __CIGAR_SOFT_CLIP) {
                svtrek_index v = bam_endpos(aln);
                if (abs_diff(v, annotated.start) <= range)
                    if (push_index(&sg.sc_start, &sg.n_sc_start, &cap_sc_start, v) < 0)
                        break;
            }

            // left (5') soft-clip; deletion END at the alignment start
            if (bam_cigar_op(cigar[0]) == __CIGAR_SOFT_CLIP) {
                svtrek_index v = aln->core.pos;
                if (abs_diff(v, annotated.end) <= range)
                    if (push_index(&sg.sc_end, &sg.n_sc_end, &cap_sc_end, v) < 0)
                        break;
            }
        }
        sam_itr_destroy(iter);
    }
    bam_destroy1(aln);

    return sg;
}

/*
 * Cluster the signatures jointly and return the refined breakpoints.
 *
 * (1) D-op tuples are clustered in 2D: two votes join the same cluster only
 *     when BOTH breakpoints agree within `consensus_interval`, so start and end
 *     stay coupled. Among clusters that reach `consensus_min_count` votes and
 *     whose centre is within `consensus_interval_range` of the annotated
 *     (pos,end) on both axes, the MAJORITY (largest) cluster wins, ties broken
 *     by closeness to the annotation -- the range gate stops it from snapping
 *     onto a stronger but unrelated neighbouring deletion. The winning cluster's
 *     start and end are then sharpened with the soft-clip signals that agree
 *     with it.
 *
 * (2) Fallback: when no D-op cluster qualifies (e.g. the deletion is longer than
 *     a read, so only split reads see it), the one-sided soft-clip signals are
 *     clustered independently around each annotated breakpoint via consensus_pos.
 *
 * Returns the supporting read count (0 = nothing qualified) and writes the
 * refined breakpoints to *ref_start / *ref_end.
 */
static size_t cluster_deletion_signatures(const del_signatures *sg, interval annotated, t_arg *params, svtrek_index *ref_start, svtrek_index *ref_end) {
    svtrek_index tol = params->consensus_interval;
    svtrek_index range = params->consensus_interval_range;
    uint32_t min_count = params->consensus_min_count;

    // majority joint cluster over D-op tuples
    uint32_t best_count = 0;
    int64_t best_dist = 0x7fffffffffffffffLL;
    svtrek_index cand_start = 0, cand_end = 0;

    *ref_start = 0;
    *ref_end = 0;

    for (size_t i = 0; i < sg->n_tuple; i++) {
        uint32_t count = 0;
        int64_t sum_start = 0, sum_end = 0;

        for (size_t j = 0; j < sg->n_tuple; j++) {
            if (abs_diff(sg->tup_start[j], sg->tup_start[i]) <= tol &&
                abs_diff(sg->tup_end[j], sg->tup_end[i])   <= tol) {
                count++;
                sum_start += sg->tup_start[j];
                sum_end   += sg->tup_end[j];
            }
        }

        if (count < min_count)
            continue;

        svtrek_index cs = (svtrek_index)((sum_start + count / 2) / count);
        svtrek_index ce = (svtrek_index)((sum_end   + count / 2) / count);

        if (abs_diff(annotated.start, cs) >= range || abs_diff(annotated.end, ce) >= range)
            continue;

        int64_t dist = (int64_t)abs_diff(annotated.start, cs) + (int64_t)abs_diff(annotated.end, ce);

        /* majority first; closeness to the annotation breaks ties. Swap the two
         * keys to prefer the closest cluster instead (consensus_pos policy). */
        if (count > best_count || (count == best_count && dist < best_dist)) {
            best_count = count;
            best_dist  = dist;
            cand_start = cs;
            cand_end   = ce;
        }
    }

    if (best_count >= min_count) {
        /* Sharpen each side with the tuple members and any soft-clip signals
         * consistent with the winning cluster. The event stays defined by the
         * joint tuple cluster; soft-clips only refine the two coordinates. */
        int64_t sum_start = 0, sum_end = 0;
        uint32_t n_start = 0, n_end = 0;

        for (uint32_t j = 0; j < sg->n_tuple; j++) {
            if (abs_diff(sg->tup_start[j], cand_start) <= tol &&
                abs_diff(sg->tup_end[j], cand_end)   <= tol) {
                sum_start += sg->tup_start[j]; n_start++;
                sum_end   += sg->tup_end[j];   n_end++;
            }
        }
        for (uint32_t j = 0; j < sg->n_sc_start; j++)
            if (abs_diff(sg->sc_start[j], cand_start) <= tol) { sum_start += sg->sc_start[j]; n_start++; }
        for (uint32_t j = 0; j < sg->n_sc_end; j++)
            if (abs_diff(sg->sc_end[j], cand_end) <= tol) { sum_end += sg->sc_end[j]; n_end++; }

        *ref_start = (svtrek_index)((sum_start + n_start / 2) / n_start);
        *ref_end   = (svtrek_index)((sum_end   + n_end   / 2) / n_end);
        return best_count;
    }

    // split-read-only fallback: cluster each side independently
    svtrek_index s = consensus_pos(sg->sc_start, sg->n_sc_start, (int)annotated.start, min_count, tol, range);
    svtrek_index e = consensus_pos(sg->sc_end,   sg->n_sc_end,   (int)annotated.end,   min_count, tol, range);
    if (s > 0 && e > 0 && e > s) {
        *ref_start = s;
        *ref_end   = e;
        size_t sc = 0, ec = 0;
        for (size_t j = 0; j < sg->n_sc_start; j++) if (abs_diff(sg->sc_start[j], s) <= tol) sc++;
        for (size_t j = 0; j < sg->n_sc_end;   j++) if (abs_diff(sg->sc_end[j], e) <= tol) ec++;
        return sc < ec ? sc : ec;
    }

    return 0;
}

/* ------------------------------------------------------------------ */
/*  Public SV entry points                                            */
/* ------------------------------------------------------------------ */

void deletion(int chrom, interval begin, interval end, interval sv_inter, t_arg *params, interval *res_inter, sv_consensus *cons) {

    res_inter->start = 0xFFFFFFFF;
    res_inter->end   = 0xFFFFFFFF;
    cons->seq = NULL; cons->len = 0; cons->support = 0;

    // collect per-read breakpoint signatures across the whole search span
    interval win = { begin.start, end.end };
    del_signatures sg = collect_deletion_signatures(chrom, win, sv_inter, params);

    // cluster them jointly and take the majority.
    svtrek_index s = 0, e = 0;
    size_t support = cluster_deletion_signatures(&sg, sv_inter, params, &s, &e);
    free_del_signatures(&sg);

    if (support > 0) {
        res_inter->start = (uint32_t)s;
        res_inter->end   = (uint32_t)e;
    }

    // base-level consensus
    svtrek_index lo = (s > 0) ? s : sv_inter.start;
    svtrek_index hi = (e > 0) ? e : sv_inter.end;
    if (hi > lo) {
        svtrek_index flank = __SV_MIN_LENGTH * 2;   /* bases of flank on each side */
        char **seqs = NULL; int *lens = NULL;
        size_t nseq = collect_spanning_seqs(chrom, lo, hi, flank, params, &seqs, &lens, 200);
        if (nseq >= params->consensus_min_count)
            run_msa(seqs, lens, nseq, cons);
        for (size_t i = 0; i < nseq; i++) free(seqs[i]);
        free(seqs); free(lens);
    }
}

void insertion(int chrom, interval begin, svtrek_index pos, t_arg *params, svtrek_index *res_start, sv_consensus *cons) {
    svtrek_index cp = refine_ins(chrom, begin, pos, params);
    *res_start = (cp == 0) ? 0xFFFFFFFF : cp;

    cons->seq = NULL; cons->len = 0; cons->support = 0;

    /* Build a consensus of the inserted (ALT) sequences. */
    if (cp > 0) {
        char **seqs = NULL; int *lens = NULL;
        size_t nseq = collect_insertion_seqs(chrom, cp, params->consensus_interval_range, params, &seqs, &lens, 500);
        if (nseq >= 1)
            run_msa(seqs, lens, nseq, cons);
        for (size_t i = 0; i < nseq; i++) free(seqs[i]);
        free(seqs); free(lens);
    }
}

void inversion(int chrom, interval begin, interval end, interval sv_inter, t_arg *params, interval *res_inter) {
    svtrek_index s = refine_point(SV_INV, chrom, begin, sv_inter.start, params);
    svtrek_index e = refine_point(SV_INV, chrom, end, sv_inter.end, params);
    res_inter->start = (s == 0) ? 0xFFFFFFFF : s;
    res_inter->end   = (e == 0) ? 0xFFFFFFFF : e;
}