#ifndef __MSA_H__
#define __MSA_H__

#include <stddef.h>

/**
 * @brief Consensus sequence produced by a multiple sequence alignment (abPOA).
 *
 * `seq` is a heap-allocated, null-terminated nucleotide string (A/C/G/T/N)
 * or NULL when no consensus could be built. `len` is its length and
 * `support` is the number of input sequences that were fed into the MSA.
 */
typedef struct _sv_consensus {
    char *seq;       /**< Consensus sequence (malloc'd, null-terminated) or NULL. */
    size_t len;      /**< Length of the consensus sequence (0 if none). */
    size_t support;  /**< Number of reads/sequences used to build the consensus. */
} sv_consensus;

/**
 * @brief Run an abPOA multiple sequence alignment and extract the consensus.
 *
 * The routine is self-contained and thread-safe: it creates and frees its own
 * abPOA handles on every call, so it may be invoked concurrently from worker
 * threads (each on independent data).
 *
 * @param seqs     Array of `n_seqs` null-terminated nucleotide strings.
 * @param lens     Array of `n_seqs` sequence lengths.
 * @param n_seqs   Number of input sequences (must be >= 1).
 * @param out      Output consensus. On success out->seq is set (caller frees
 *                 it with sv_consensus_free); on failure out->seq stays NULL.
 *
 * @return 0 on success (a consensus was produced), -1 otherwise.
 */
int run_msa(char **seqs, int *lens, size_t n_seqs, sv_consensus *out);

/**
 * @brief Release the memory held by an sv_consensus produced by run_msa.
 */
void sv_consensus_free(sv_consensus *c);

#endif
