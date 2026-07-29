#include "msa.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <pthread.h>
#include <abpoa.h>

// Nucleotide -> 2-bit code table: A/a=0 C/c=1 G/g=2 T/t=3, everything else=4 (N).

static unsigned char nt4_table[256];
static pthread_once_t nt4_once = PTHREAD_ONCE_INIT;

static void init_nt4(void) {
    for (int i = 0; i < 256; i++) nt4_table[i] = 4;
    nt4_table[(unsigned char)'A'] = nt4_table[(unsigned char)'a'] = 0;
    nt4_table[(unsigned char)'C'] = nt4_table[(unsigned char)'c'] = 1;
    nt4_table[(unsigned char)'G'] = nt4_table[(unsigned char)'g'] = 2;
    nt4_table[(unsigned char)'T'] = nt4_table[(unsigned char)'t'] = 3;
}

// code -> base for the consensus (0..3 => ACGT, anything else => N).
static const char code2base[5] = { 'A', 'C', 'G', 'T', 'N' };

int run_msa(char **seqs, int *lens, int n_seqs, sv_consensus *out) {
    out->seq = NULL;
    out->len = 0;
    out->support = 0;

    if (n_seqs < 1)
        return -1;

    pthread_once(&nt4_once, init_nt4);

    abpoa_t *ab = abpoa_init();
    abpoa_para_t *abpt = abpoa_init_para();

    /* We only need the consensus, not the row-column MSA dump. */
    abpt->out_msa = 0;
    abpt->out_cons = 1;
    abpt->max_n_cons = 1;
    abpt->disable_seeding = 1;   /* SV sequences are short; use full POA, not minimizer partitioning. */
    abpt->progressive_poa = 0;
    abpt->verbose = 0;
    abpoa_post_set_para(abpt);

    /* Encode the input sequences to abPOA's 2-bit representation. */
    uint8_t **bseqs = (uint8_t **)malloc(sizeof(uint8_t *) * n_seqs);
    if (bseqs == NULL) {
        abpoa_free(ab);
        abpoa_free_para(abpt);
        return -1;
    }
    for (int i = 0; i < n_seqs; i++) {
        bseqs[i] = (uint8_t *)malloc(lens[i] > 0 ? lens[i] : 1);
        for (int j = 0; j < lens[i]; j++)
            bseqs[i][j] = nt4_table[(unsigned char)seqs[i][j]];
    }

    // No per-base quality weights (NULL), no output stream (NULL).
    abpoa_msa(ab, abpt, n_seqs, NULL, lens, bseqs, NULL, NULL);

    abpoa_cons_t *abc = ab->abc;
    int rc = -1;
    if (abc != NULL && abc->n_cons > 0 && abc->cons_len[0] > 0) {
        int L = abc->cons_len[0];
        char *s = (char *)malloc(L + 1);
        if (s != NULL) {
            for (int j = 0; j < L; j++) {
                uint8_t b = abc->cons_base[0][j];
                s[j] = (b < 4) ? code2base[b] : 'N';
            }
            s[L] = '\0';
            out->seq = s;
            out->len = L;
            out->support = n_seqs;
            rc = 0;
        }
    }

    for (int i = 0; i < n_seqs; i++) free(bseqs[i]);
    free(bseqs);
    abpoa_free(ab);
    abpoa_free_para(abpt);

    return rc;
}

void sv_consensus_free(sv_consensus *c) {
    if (c != NULL && c->seq != NULL) {
        free(c->seq);
        c->seq = NULL;
        c->len = 0;
        c->support = 0;
    }
}
