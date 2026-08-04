#include "refinement.h"
#include "msa.h"
#include <htslib/sam.h>


#define abs_int(a) (a < 0 ? -(a) : a)


/* ------------------------------------------------------------------ */
/*  Consensus of a set of candidate breakpoint locations              */
/* ------------------------------------------------------------------ */

/*
 * Legacy, currently-unused helper kept for reference. Finds the value that
 * anchors the largest cluster (values within `consensus_interval`).
 */
int consensus(int *arr, int size, int consensus_min_count, int consensus_interval) {
    int consensus_val = -1;
    int max_count = consensus_min_count - 1;

    quicksort(arr, 0, size-1);

    for (int i = 0; i < size; i++) {
        int count = 1;
        // count how many values fall within consensus_interval of arr[i]
        for (int j = i + 1; j < size && arr[j] <= arr[i] + consensus_interval; j++) {
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
int consensus_pos(int *locations, int size, int pos, int consensus_min_count, int consensus_interval, int consensus_interval_range) {

    if (size < consensus_min_count) {
        return -1;
    }

    quicksort(locations, 0, size - 1);

    int best_candidate = -1;
    int best_distance = 0x7fffffff;
    int best_count = 0;

    for (int i = 0; i < size; i++) {
        /* Grow a cluster of locations within consensus_interval of one
         * another, anchored at locations[i] (array is sorted). */
        int count = 1;
        int64_t total = locations[i];
        for (int j = i + 1; j < size && locations[j] <= locations[i] + consensus_interval; j++) {
            count++;
            total += locations[j];
        }

        if (count < consensus_min_count)
            continue;

        int candidate = (int)((total + count / 2) / count);
        int distance = abs_int(pos - candidate);

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

int refine_start(sv_type_t sv_type, int chrom, interval inter, uint32_t imprecise_pos, t_arg *params) {

    int capacity = 100;
    int size = 0;
    int *start_locations = (int *)malloc(sizeof(int)*capacity);
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
                        int *temp = (int *)realloc(start_locations, sizeof(int)*capacity);
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
                    int *temp = (int *)realloc(start_locations, sizeof(int)*capacity);
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
    int result = consensus_pos(start_locations, size, imprecise_pos, params->consensus_min_count, params->consensus_interval, params->consensus_interval_range);
    free(start_locations);
    return result;
}

int refine_end(sv_type_t sv_type, int chrom, interval inter, uint32_t imprecise_pos, t_arg *params) {

    int capacity = 100;
    int size = 0;
    int *end_locations = (int *)malloc(sizeof(int)*capacity);
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
                        int *temp = (int *)realloc(end_locations, sizeof(int)*capacity);
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
                    int *temp = (int *)realloc(end_locations, sizeof(int)*capacity);
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
    int result = consensus_pos(end_locations, size, imprecise_pos, params->consensus_min_count, params->consensus_interval, params->consensus_interval_range);
    free(end_locations);
    return result;
}

int refine_point(sv_type_t sv_type, int chrom, interval inter, uint32_t imprecise_pos, t_arg *params) {

    int capacity = 100;
    int size = 0;
    int *locations = (int *)malloc(sizeof(int)*capacity);
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
                        int *temp = (int *)realloc(locations, sizeof(int)*capacity);
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
    int result = consensus_pos(locations, size, imprecise_pos, params->consensus_min_count, params->consensus_interval, params->consensus_interval_range);
    free(locations);
    return result;
}

int refine_ins(int chrom, interval inter, uint32_t imprecise_pos, t_arg *params) {

    int capacity = 100;
    int size = 0;
    int *locations = (int *)malloc(sizeof(int)*capacity);
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
                        int *temp = (int *)realloc(locations, sizeof(int)*capacity);
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
    int result = consensus_pos(locations, size, imprecise_pos, params->consensus_min_count, params->consensus_interval, params->consensus_interval_range);
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
static int collect_insertion_seqs(int chrom, int cons_pos, int tol, t_arg *params,
                                  char ***seqs_out, int **lens_out, int max_seqs) {
    char **seqs = (char **)malloc(sizeof(char *) * max_seqs);
    int   *lens = (int  *)malloc(sizeof(int)    * max_seqs);
    int n = 0;

    if (seqs == NULL || lens == NULL) {
        free(seqs); free(lens);
        *seqs_out = NULL; *lens_out = NULL;
        return 0;
    }

    int lo = cons_pos - tol; if (lo < 1) lo = 1;
    int hi = cons_pos + tol;

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
                int len = bam_cigar_oplen(cigar[i]);

                if (op == __CIGAR_INSERTION && len >= __SV_MIN_LENGTH && abs_int((int)ref - cons_pos) <= tol) {
                    char *s = (char *)malloc(len + 1);
                    if (s != NULL) {
                        for (int j = 0; j < len; j++)
                            s[j] = seq_nt16_str[bam_seqi(qseq, qpos + j)];
                        s[len] = '\0';
                        seqs[n] = s;
                        lens[n] = len;
                        n++;
                    }
                }

                if (bam_cigar_type(op) & 1) qpos += len;   /* consumes query */
                if (bam_cigar_type(op) & 2) ref  += len;   /* consumes reference */

                if ((int)ref > hi + len) break;
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
static int collect_spanning_seqs(int chrom, int lo, int hi, int flank, t_arg *params,
                                 char ***seqs_out, int **lens_out, int max_seqs) {
    int win_lo = lo - flank; if (win_lo < 0) win_lo = 0;
    int win_hi = hi + flank;

    char **seqs = (char **)malloc(sizeof(char *) * max_seqs);
    int   *lens = (int  *)malloc(sizeof(int)    * max_seqs);
    int n = 0;

    if (seqs == NULL || lens == NULL) {
        free(seqs); free(lens);
        *seqs_out = NULL; *lens_out = NULL;
        return 0;
    }

    int cap_seq = (win_hi - win_lo) + 16;
    int tid = chrom;
    bam1_t *aln = bam_init1();
    hts_itr_t *iter = sam_itr_queryi(params->hargs.bam_file_index, tid, win_lo, win_hi);
    if (iter) {
        while (n < max_seqs && sam_itr_next(params->hargs.fp_in, iter, aln) > 0) {
            int aln_start = aln->core.pos;
            int aln_end   = bam_endpos(aln);           /* 0-based, exclusive */
            /* require the read to fully span the window */
            if (aln_start > win_lo || aln_end < win_hi)
                continue;

            uint32_t ref = aln->core.pos;
            uint32_t qpos = 0;
            uint32_t *cigar = bam_get_cigar(aln);
            uint8_t  *qseq  = bam_get_seq(aln);

            char *s = (char *)malloc(cap_seq + 1);
            if (s == NULL) continue;
            int sl = 0;

            for (uint32_t i = 0; i < aln->core.n_cigar; i++) {
                int op  = bam_cigar_op(cigar[i]);
                int len = bam_cigar_oplen(cigar[i]);
                int cq  = bam_cigar_type(op) & 1;      /* consumes query */
                int cr  = bam_cigar_type(op) & 2;      /* consumes reference */

                if (cq && cr) {                        /* M / = / X */
                    for (int j = 0; j < len && sl < cap_seq; j++) {
                        int rp = (int)ref + j;
                        if (rp >= win_lo && rp < win_hi)
                            s[sl++] = seq_nt16_str[bam_seqi(qseq, qpos + j)];
                    }
                } else if (op == __CIGAR_INSERTION) {  /* inserted bases anchored at ref */
                    if ((int)ref >= win_lo && (int)ref < win_hi)
                        for (int j = 0; j < len && sl < cap_seq; j++)
                            s[sl++] = seq_nt16_str[bam_seqi(qseq, qpos + j)];
                }

                if (cq) qpos += len;
                if (cr) ref  += len;
            }

            if (sl > 0) {
                s[sl] = '\0';
                seqs[n] = s;
                lens[n] = sl;
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
/*  Public SV entry points                                            */
/* ------------------------------------------------------------------ */

void deletion(int chrom, interval begin, interval end, interval sv_inter, t_arg *params, interval *res_inter, sv_consensus *cons) {
    int s = refine_start(SV_DEL, chrom, begin, sv_inter.start, params);
    int e = refine_end(SV_DEL, chrom, end, sv_inter.end, params);

    res_inter->start = (s < 0) ? 0xFFFFFFFF : (uint32_t)s;
    res_inter->end   = (e < 0) ? 0xFFFFFFFF : (uint32_t)e;

    cons->seq = NULL; cons->len = 0; cons->support = 0;

    /* Build a consensus of the reads spanning the deletion. */
    int lo = (s >= 0) ? s : (int)sv_inter.start;
    int hi = (e >= 0) ? e : (int)sv_inter.end;
    if (hi > lo) {
        int flank = __SV_MIN_LENGTH * 2;   /* bases of flank on each side */
        char **seqs = NULL; int *lens = NULL;
        int nseq = collect_spanning_seqs(chrom, lo, hi, flank, params, &seqs, &lens, 200);
        if (nseq >= params->consensus_min_count)
            run_msa(seqs, lens, nseq, cons);
        for (int i = 0; i < nseq; i++) free(seqs[i]);
        free(seqs); free(lens);
    }
}

void insertion(int chrom, interval begin, uint32_t pos, t_arg *params, uint32_t *res_start, sv_consensus *cons) {
    int cp = refine_ins(chrom, begin, pos, params);
    *res_start = (cp < 0) ? 0xFFFFFFFF : (uint32_t)cp;

    cons->seq = NULL; cons->len = 0; cons->support = 0;

    /* Build a consensus of the inserted (ALT) sequences. */
    if (cp >= 0) {
        char **seqs = NULL; int *lens = NULL;
        int nseq = collect_insertion_seqs(chrom, cp, params->consensus_interval_range, params, &seqs, &lens, 500);
        if (nseq >= 1)
            run_msa(seqs, lens, nseq, cons);
        for (int i = 0; i < nseq; i++) free(seqs[i]);
        free(seqs); free(lens);
    }
}

void inversion(int chrom, interval begin, interval end, interval sv_inter, t_arg *params, interval *res_inter) {
    int s = refine_point(SV_INV, chrom, begin, sv_inter.start, params);
    int e = refine_point(SV_INV, chrom, end, sv_inter.end, params);
    res_inter->start = (s < 0) ? 0xFFFFFFFF : (uint32_t)s;
    res_inter->end   = (e < 0) ? 0xFFFFFFFF : (uint32_t)e;
}