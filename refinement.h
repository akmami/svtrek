#ifndef __REFINEMENT_H__
#define __REFINEMENT_H__

#include "params.h"
#include "msa.h"
#include <stdlib.h>

#define abs(a) (a < 0 ? -(a) : a)

/**
 * @brief Processes a deletion event within a given chromosomal region.
 *
 * This function refines and validates the presence of a structural variant (SV) deletion 
 * within the specified chromosomal interval. If no refinement is performed, the function 
 * returns -1. Otherwise, it updates the result interval with refined breakpoints.
 *
 * @param chrom The chromosome identifier.
 * @param begin The starting interval of the deletion.
 * @param end The ending interval of the deletion.
 * @param sv_inter The structural variant interval to be refined.
 * @param params Additional parameters for processing the deletion.
 * @param res_inter Pointer to an interval that stores the refined deletion boundaries.
 * @param cons Output consensus of the reads spanning the deletion (abPOA MSA).
 *
 * @return -1 if no refinement is done; otherwise, res_inter is updated with refined coordinates.
 */
void deletion(int chrom, interval begin, interval end, interval sv_inter, t_arg *params, interval *res_inter, sv_consensus *cons);

/**
 * @brief Processes an insertion event within a given chromosomal region.
 *
 * This function refines and validates the presence of a structural variant (SV) insertion 
 * within the specified chromosomal position. If no refinement is performed, the function 
 * returns -1. Otherwise, it updates the result position with a refined breakpoint.
 *
 * @param chrom The chromosome identifier.
 * @param begin The interval containing the insertion position.
 * @param pos The specific position where the insertion occurs.
 * @param params Additional parameters for processing the insertion.
 * @param res_start Pointer to a variable that stores the refined insertion position.
 * @param cons Output consensus of the inserted (ALT) sequences (abPOA MSA).
 *
 * @return -1 if no refinement is done; otherwise, res_start is updated with the refined position.
 */
void insertion(int chrom, interval begin, uint32_t pos, t_arg *params, uint32_t *res_start, sv_consensus *cons);

/**
 * @brief Processes an inversion event within a given chromosomal region.
 *
 * This function refines and validates the presence of a structural variant (SV) inversion 
 * within the specified chromosomal interval. If no refinement is performed, the function 
 * returns -1. Otherwise, it updates the result interval with refined breakpoints.
 *
 * @param chrom The chromosome identifier.
 * @param begin The starting interval of the inversion.
 * @param end The ending interval of the inversion.
 * @param sv_inter The structural variant interval to be refined.
 * @param params Additional parameters for processing the inversion.
 * @param res_inter Pointer to an interval that stores the refined inversion boundaries.
 *
 * @return -1 if no refinement is done; otherwise, res_inter is updated with refined coordinates.
 */
void inversion(int chrom, interval begin, interval end, interval sv_inter, t_arg *params, interval *res_inter);

void quicksort(int *array, int low, int high);

#endif
