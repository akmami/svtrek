#include <stdlib.h>
#include <stdio.h>
#include <htslib/sam.h>
#include "sliding_window.h"
#include "refinement.h"


int refine_ins_disc(int chrom, interval inter, t_arg *params, int windowSize, int slideSize) {
    int bestCandidateOverall = -1; 
    int maxSupportOverall = 0;

    for (uint32_t sub_start = inter.start; sub_start < inter.end; sub_start += windowSize) {
        uint32_t sub_end = sub_start + windowSize;
        if (sub_end > inter.end)
            sub_end = inter.end;  

        
        int capacity = 100; // dynamic but initialized with 100
        int size = 0;
        int *locations = (int *)malloc(sizeof(int) * capacity);
        if (locations == NULL) {
            fprintf(stderr, "[ERROR] Couldn't allocate array for positions in sub-window [%d, %d].\n", sub_start, sub_end);
            return -1;
        }

        bam1_t *aln = bam_init1();
        hts_itr_t *iter = sam_itr_queryi(params->hargs.bam_file_index, chrom - 1, sub_start - 1, sub_end - 1);
        if (iter) {
            while (sam_itr_next(params->hargs.fp_in, iter, aln) > 0) {
                uint32_t reference_pos = aln->core.pos;
                uint32_t *cigar = bam_get_cigar(aln);
                for (uint32_t i = 0; i < aln->core.n_cigar; i++) {
                    if (bam_cigar_op(cigar[i]) == __CIGAR_INSERTION &&
                        bam_cigar_oplen(cigar[i]) >= __SV_MIN_LENGTH) {
                        if (size == capacity) {
                            capacity = (int)(capacity * 1.5);
                            int *temp = (int *)realloc(locations, sizeof(int) * capacity);
                            if (temp == NULL) {
                                fprintf(stderr, "[ERROR] Couldn't reallocate locations array.\n");
                                free(locations);
                                bam_destroy1(aln);
                                return -1;
                            }
                            locations = temp;
                        }
                        locations[size++] = reference_pos;
                    }
                    if (bam_cigar_op(cigar[i]) != __CIGAR_INSERTION &&
                        bam_cigar_op(cigar[i]) != __CIGAR_SOFT_CLIP) {
                        reference_pos += bam_cigar_oplen(cigar[i]);
                    }
                    if (reference_pos > sub_end)
                        break;
                }
            }
            sam_itr_destroy(iter);
        }
        bam_destroy1(aln);

        if (size == 0) {
            free(locations);
            continue;
        }

        quicksort(locations, 0, size - 1);

        int bestCandidate = -1;
        int maxSupport = 0;
        
        for (int i = 0; i < size; i += slideSize) {
            int end = i;
            while (end < size && (locations[end] - locations[i]) <= windowSize) {
                end++;
            }
            int support = end - i; 
            if (support >= params->consensus_min_count && support > maxSupport) {
                maxSupport = support;
                int sum = 0;
                for (int j = i; j < end; j++) {
                    sum += locations[j];
                }
                bestCandidate = (sum + support / 2) / support; 
            }
        }

        if (bestCandidate != -1) {
            printf("INS Discovery in window [%d, %d] at position %d with support %d\n", sub_start, sub_end, bestCandidate, maxSupport);
            if (maxSupport > maxSupportOverall) {
                maxSupportOverall = maxSupport;
                bestCandidateOverall = bestCandidate;
            }
        }
        free(locations);
    }

    return bestCandidateOverall;
}


