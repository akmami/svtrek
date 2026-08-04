#ifndef __UTILS_H__
#define __UTILS_H__

#include "params.h"
#include <stdint.h>
#include <string.h>


#define validate(temp) \
		if (!temp) { \
			fprintf(stderr, "[ERROR] Memory allocation failed\n"); \
			exit(1); \
		}

#define available(type, ptr, size, capacity) \
		if (size == capacity) { \
			capacity = capacity * 2; \
			type *temp = realloc(ptr, capacity * sizeof(type)); \
			validate(temp) \
			ptr = temp; \
		}

/**
 * @brief Quick sort for int.
 * 
 * @param array The array to be sorted
 * @param low   The left most element's index
 * @param high  The right most element's index
 */
void quicksort(int *array, int low, int high);

/**
 * @brief Reverse the given array.
 * 
 * @param arr   The array to be reverse.
 * @param size  The size of the array.
 */
void reverse(uint64_t *arr, int size);

/**
 * @brief Correct the start and end indices of the given elements.abort
 * 
 * In some GAF files, the start and end indices of the reverse complement
 * reads are not corrected w.r.t backbone reference genome. Hence, the start
 * and end indices needs to be swapped.
 * 
 * @param start The start index
 * @param end   The end index
 * @param len   The lenght of the entire region
 * @param new_start The new start based on the correction
 * @param new_end   The new end based on the correction
 */
void fix_indices(int start, int end, int len, int *new_start, int *new_end);

/**
 * @brief Parses a CIGAR string into a linear sequence of operations.
 *
 * This function expands the CIGAR string (e.g., "5M2I3D") into an array of individual operations 
 * (e.g., ['M','M','M','M','M','I','I','D','D','D']).
 *
 * @param cigar     The input CIGAR string (null-terminated).
 * @param ops       Output buffer to hold expanded operations.
 * @param max_ops   Maximum number of operations `ops` can hold.
 * @param rev       If non-zero, the output ops array will be reversed (useful for reverse strand).
 *
 * @return The number of parsed operations on success, 
 *         or -1 if the number of operations exceeds `max_ops`.
 */
int parse_cigar(char *cigar, char *ops, int max_ops, int rev);

/**
 * @brief Parses a node path string and extracts the strand and node ID.
 *
 * This function reads the strand character ('>' or '<') followed by the numeric ID from a VG-style path string 
 * (e.g., ">1234", "<5678").
 *
 * @param path      Pointer to the path string
 * @param id        Output pointer to store the parsed node ID.
 * @param strand    Output pointer to store the strand character ('>' or '<').
 *
 * @return The number of characters consumed from the path string,
 *         or 0 if the path is empty.
 */
int next_node(const char *path, uint64_t *id, char *strand);

FILE *io_open(const char* file_path, char **line, int cap, const char *mode);
void io_close(FILE *file, char **str);
int io_read(FILE *file, char **str, size_t *cap);
void free_segments(segment **segments, int segment_size);

#endif