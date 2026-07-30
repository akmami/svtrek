#ifndef PARAMS_H
#define PARAMS_H

#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <sys/stat.h>
#include <pthread.h>

//                                Consumes     query   reference   Op
#define __CIGAR_ALIGNMENT_MATCH       0     //  yes     yes         M
#define __CIGAR_INSERTION             1     //  yes     no          I
#define __CIGAR_DELETION              2     //  no      yes         D
#define __CIGAR_SKIPPED               3     //  no      yes         N
#define __CIGAR_SOFT_CLIP             4     //  yes     no          S
#define __CIGAR_HARD_CLIP             5     //  no      no          H
#define __CIGAR_PADDING               6     //  no      no          P
#define __CIGAR_SEQUENCE_MATCH        7     //  yes     yes         =
#define __CIGAR_SEQUENCE_MISMATCH     8     //  yes     yes         X

#define CIGAR_QUE(x) (x==__CIGAR_ALIGNMENT_MATCH || x==__CIGAR_INSERTION || x==__CIGAR_SOFT_CLIP || x==__CIGAR_SEQUENCE_MATCH || x==__CIGAR_SEQUENCE_MISMATCH )
#define CIGAR_REF(x) (x==__CIGAR_ALIGNMENT_MATCH || x==__CIGAR_DELETION || x==__CIGAR_SKIPPED || x==__CIGAR_SEQUENCE_MATCH || x==__CIGAR_SEQUENCE_MISMATCH )

#define __FLAG_MULTIPLE_SEGMENTS          0x1
#define __FLAG_SECONDARY_ALIGNMENT        0x100
#define __FLAG_SUPPLEMENTARY_ALIGNMENT    0x800

#define __WIDER_INTERVAL           20000
#define __MEDIAN_INTERVAL          10000
#define __NARROW_INTERVAL          2000
#define __CONSENSUS_INTERVAL_RANGE 150   // max bp a breakpoint may be relocated in refinement (CLI: --consensus-interval-range)
#define __CONSENSUS_INTERVAL       5
#define __CONSENSUS_MIN_COUNT      3
#define __SV_MIN_LENGTH            50
//#define CONSENSUS_COUNT_PERC     0.3

#define __THREAD_NUMBER           4
#define __THREAD_POOL_LOAD_FACTOR 2

#define MAX_LINE 1048576
#define MAX_CIGAR 131072
#define MAX_SEQ 1048576


typedef enum _mode {
    MODE_DISC,
    MODE_AUDT
} program_mode;

typedef struct _audt_args {
    // input arguments
    const char* bam_file;
    const char* vcf_file;
    const char* output_file;
    int thread_number;
    int verbose;
    int tload_factor;
    // program arguments
    int wider_interval;
    int median_interval;
    int narrow_interval;
    int consensus_interval_range;
    int consensus_interval;
    int consensus_min_count;
    program_mode mode;
} audt_args;

typedef struct _line_queue {
    char **lines;   /**< Queue to store lines extracted from VCF file for threads. */
    int size;       /**< The size of the queue. */
    int capacity;   /**< Capacity of the queue. */
    int front;      /**< The index for the pushing point. */
    int rear;       /**< The index for the popping point. */
} line_queue;

typedef struct _htslib_params {
    samFile *fp_in;
    bam_hdr_t *bam_hdr;
    hts_idx_t *bam_file_index;
} htslib_params;

typedef struct _t_arg {
    int wider_interval;
    int median_interval;
    int narrow_interval;
    int consensus_interval_range;
    int consensus_interval;
    int consensus_min_count;
    htslib_params hargs;
    line_queue *queue;
    pthread_mutex_t *queue_mutex;
    pthread_mutex_t *out_err_mutex;
    pthread_cond_t *cond_not_full;
    pthread_cond_t *cond_not_empty;
    int *exit_signal;
} t_arg;

typedef struct _disc_args {
    // input arguments
    const char* gfa_file;
    const char* gaf_file;
    const char* fq_file;
    const char* output_file;
    int thread_number;
    int verbose;
    int tload_factor;
    // program arguments
    int consensus_interval_range;
    int consensus_interval;
    int consensus_min_count;
    program_mode mode;
} disc_args;

typedef enum {
    SV_UNKNOWN = 0,
    SV_INS,    // Insertion
    SV_DEL,    // Deletion
    SV_INV,    // Inversion
    SV_DUP,    // Duplication
    SV_TRA,    // Translocation
    SV_BND     // Breakend
} sv_type_t;

typedef struct _interval {
    uint32_t start;
    uint32_t end;
} interval;

typedef struct _alignment {
    char *readName;
    int readLen;
    int readStart;      // The start position of the query (read) in the alignment.
    int readEnd;        // The end position of the query (read) in the alignment.
    char strand;        // + or -
    char *path;         // node names in < or > directions
    int pathLen;
    int pathStart;      // The start position of the alignment on the reference node path.
    int pathEnd;        // The end position of the alignment on the reference node path.
    int matches;        // The number of matching bases between the query and the reference.
    int blockLen;       // The length of the aligned block (number of bases in the alignment block).
    int qual;           // The mapping quality score, 255 not avaialble 
    int xdi;            // NM:i:22 The number of mismatches, deletions, and insertions in the alignment (represented as an integer value).
    float score;        // AS:f:14239.3 (optional) The alignment score as a floating-point number, if provided (this field is optional).
    float divergence;   // dv:f:0.00153685 The divergence score, which quantifies how different the alignment is by calculating a ratio of mismatches, deletions, and insertions relative to the total alignment length.
    float identity;     // id:f:0.998463 The identity score, which is the ratio of matches to the total alignment length (matches, mismatches, deletions, and insertions).
    char *cigar;        // cg:Z:3=1I (optional)
} alignment;

typedef struct _gaf_breakpoint {
    int readStart;
    int readEnd;
    int rc;
    int offset;
    sv_type_t type;
} gaf_breakpoint;

typedef struct segment {
    uint64_t id;
    char *seq;
    int start;
    int end;
    int rank;
    struct segment *next;
} segment;

#endif // PARAMS_H
