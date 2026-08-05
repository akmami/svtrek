#include "audit.h"

static inline sv_type_t parse_sv_type(const char *sv_str) {
    if (strcmp(sv_str, "INS") == 0 || strcmp(sv_str, "INS:ME") == 0 ) return SV_INS;
    if (strcmp(sv_str, "DEL") == 0 || strcmp(sv_str, "DEL:ME") == 0 ) return SV_DEL;
    if (strcmp(sv_str, "INV") == 0) return SV_INV;
    if (strcmp(sv_str, "DUP") == 0) return SV_DUP;
    if (strcmp(sv_str, "TRA") == 0) return SV_TRA;
    if (strcmp(sv_str, "BND") == 0) return SV_BND;
    return SV_UNKNOWN;
}

static inline char * sv_to_str(sv_type_t type) {
    if (type == SV_INS) return "SV_INS";
    else if (type == SV_DEL) return "SV_DEL";
    else if (type == SV_INV) return "SV_INV";
    else if (type == SV_DUP) return "SV_DUP";
    else if (type == SV_TRA) return "SV_TRA";
    else if (type == SV_BND) return "SV_BND";
    else return "UNKNOWN";
}

void line_queue_init(line_queue *queue, int capacity) {
    queue->lines = (char **)malloc(capacity * sizeof(char *));
    queue->size = 0;
    queue->capacity = capacity;
    queue->front = 0;
    queue->rear = 0;
}

static inline void line_queue_push(line_queue *queue, char *line, pthread_mutex_t *queue_mutex, pthread_cond_t *cond_not_full, pthread_cond_t *cond_not_empty) {
    pthread_mutex_lock(queue_mutex);
    while (queue->size == queue->capacity) {
        pthread_cond_wait(cond_not_full, queue_mutex);
    }
    queue->lines[queue->rear] = line;
    queue->rear = (queue->rear + 1) % queue->capacity;
    queue->size++;
    pthread_cond_signal(cond_not_empty);
    pthread_mutex_unlock(queue_mutex);
}

static inline char *line_queue_pop(line_queue *queue, pthread_mutex_t *queue_mutex, pthread_cond_t *cond_not_full, pthread_cond_t *cond_not_empty, int *exit_signal) {
    pthread_mutex_lock(queue_mutex);
    while (queue->size == 0 && *(exit_signal) == 0) {
        pthread_cond_wait(cond_not_empty, queue_mutex);
    }
    // Only stop once the queue has been fully drained. Checking exit_signal
    // alone would discard lines that were still queued when the producer
    // finished (the producer can set exit_signal before any worker pops).
    if (queue->size == 0 && *(exit_signal) == 1) {
        pthread_mutex_unlock(queue_mutex);
        return NULL;
    }
    char *line = queue->lines[queue->front];
    queue->front = (queue->front + 1) % queue->capacity;
    queue->size--;
    pthread_cond_signal(cond_not_full);
    pthread_mutex_unlock(queue_mutex);
    return line;
}

// Copy the value of an INFO sub-field (e.g. "END=" or "SVTYPE=") into dst.
// Returns 1 if the key was found, 0 otherwise. dst is always NUL-terminated,
// and truncated to fit when the value is longer than dst_size - 1.
static inline int info_field_value(const char *info, const char *key, char *dst, size_t dst_size) {
    char *value = strstr(info, key);
    if (value == NULL) {
        return 0;
    }
    value += strlen(key);

    char *value_end = strchr(value, ';');
    size_t value_len = value_end ? (size_t)(value_end - value) : strlen(value);
    if (value_len >= dst_size) {
        value_len = dst_size - 1;
    }

    memcpy(dst, value, value_len);
    dst[value_len] = '\0';
    return 1;
}

void thread_func(void *_params) {

    t_arg *targs = (t_arg *)_params;
    line_queue *queue = targs->queue;

    while (1) {
        char *line = line_queue_pop(queue, targs->queue_mutex, targs->cond_not_full, targs->cond_not_empty, targs->exit_signal);
        if (line == NULL) {
            break;
        }

        char *saveptr;

        // column 1: CHROM; resolve to a BAM target id
        char *contig = strtok_r(line, "\t", &saveptr);
        int contig_index = sam_hdr_name2tid(targs->hargs.bam_hdr, contig);
        if (contig_index < 0) {
            fprintf(stderr, "[ERROR] VCF: contig '%s' not found in BAM header, skipping\n", contig);
            free(line);
            continue;
        }

        // column 2: POS
        char *pos_str = strtok_r(NULL, "\t", &saveptr);
        if (pos_str == NULL) {
            fprintf(stderr, "[ERROR] VCF: no index at line: %s\n", line);
            free(line);
            continue;
        }

        svtrek_index pos = strtol(pos_str, NULL, 10);
        if (pos == 0 && pos_str[0] != '0') { // check for conversion error
            fprintf(stderr, "[ERROR] Variant %s, begin: %s filtered out because of PARSE-POS error.\n", contig, pos_str);
            free(line);
            continue;
        }

        strtok_r(NULL, "\t", &saveptr);         // column 3: ID (skip)

        // column 4: REF
        char *ref = strtok_r(NULL, "\t", &saveptr);
        svtrek_index ref_len = strlen(ref);

        // column 5: ALT (comma-separated); track shortest and longest allele
        char *alt = strtok_r(NULL, "\t", &saveptr);
        char *alt_saveptr;
        char *alt_token = strtok_r(alt, ",", &alt_saveptr);
        size_t max_alt_len = 0;
        size_t min_alt_len = 0x7FFFFFFF;
        while (alt_token != NULL) {
            size_t alt_len = strlen(alt_token);
            if (alt_len > max_alt_len) {
                max_alt_len = alt_len;
            }
            if (alt_len < min_alt_len) {
                min_alt_len = alt_len;
            }
            alt_token = strtok_r(NULL, ",", &alt_saveptr);
        }

        strtok_r(NULL, "\t", &saveptr);              // column 6: QUAL (skip)
        strtok_r(NULL, "\t", &saveptr);              // column 7: FILTER (skip)
        char *info = strtok_r(NULL, "\t", &saveptr); // column 8: INFO

        // END position: read from INFO "END=", otherwise POS + REF length
        svtrek_index end;
        char end_buf[32];
        if (info_field_value(info, "END=", end_buf, sizeof(end_buf))) {
            end = strtol(end_buf, NULL, 10);
            if (end == 0 && end_buf[0] != '0') {
#ifdef DEBUG
                fprintf(stderr, "[FAIL] Contig %s, begin: %u, end: %u filtered out because of PARSE-END error.\n", contig, pos, end);
#endif
                free(line);
                continue;
            }
        } else {
            end = pos + ref_len;
        }

        // SV type: read from INFO "SVTYPE=", otherwise infer from REF/ALT lengths
        sv_type_t sv_type = SV_UNKNOWN;
        char sv_buf[16];
        if (info_field_value(info, "SVTYPE=", sv_buf, sizeof(sv_buf))) {
            sv_type = parse_sv_type(sv_buf);
        } else if (ref_len == 1 && __SV_MIN_LENGTH < max_alt_len) {
            sv_type = SV_INS;
        } else if (__SV_MIN_LENGTH < ref_len && min_alt_len == 1) {
            sv_type = SV_DEL;
        }

        if (info_field_value(info, "SVLEN=", sv_buf, sizeof(sv_buf))) {
            int sv_len = atoi(sv_buf);
            if (sv_len < __SV_MIN_LENGTH) {
#ifdef DEBUG
                fprintf(stderr, "[FAIL] Contig %s, begin: %u, end: %u, len: %d filtered out because of SV len.\n", contig, pos, end, sv_len);
#endif
                free(line);
                continue;
            }
        }

        // Only INS/DEL/INV are implemented so far
        if (sv_type != SV_INS && sv_type != SV_DEL) { // && sv_type != SV_INV) {
#ifdef DEBUG
            if (sv_type != SV_UNKNOWN) {
                fprintf(stderr, "[FAIL] Skipping unsupported SV type %s at %s:%s\n", sv_to_str(sv_type), contig, pos_str);
            }
#endif
            free(line);
            continue;
        }

#ifdef DEBUG
        fprintf(stderr, "[INFO] Processing %s, begin: %u, end: %u, %s\n", contig, pos, end, sv_to_str(sv_type));
#endif

        switch (sv_type) {
        case SV_INS:
            {
                svtrek_index begin_start = pos > targs->median_interval ? pos - targs->median_interval : 1;
                interval begin = {begin_start, pos + targs->median_interval};
                svtrek_index result;
                sv_consensus cons;
                insertion(contig_index, begin, pos, targs, &result, &cons);

                pthread_mutex_lock(targs->out_err_mutex);
                if (result == 0xFFFFFFFF) {
                    printf("(INS) contig: %d, org pos: %u, ref pos: NA", contig_index, pos);
                } else {
                    printf("(INS) contig: %d, org pos: %u, ref pos: %u, diff: %d", contig_index, pos, result, result - pos);
                }
                if (cons.seq != NULL) {
                    printf(", cons_len: %ld, support: %ld, cons: %.60s%s\n", cons.len, cons.support, cons.seq, cons.len > 60 ? "..." : "");
                } else {
                    printf(", cons_len: NA\n");
                }
                pthread_mutex_unlock(targs->out_err_mutex);

                sv_consensus_free(&cons);
            }
            break;
        case SV_DEL:
            {
                svtrek_index del_begin_start = pos > targs->wider_interval ? pos - targs->wider_interval : 1;
                svtrek_index del_end_start = end > targs->narrow_interval ? end - targs->narrow_interval : 1;
                interval del_begin = {del_begin_start, pos + targs->narrow_interval};
                interval del_end = {del_end_start, end + targs->narrow_interval};
                interval sv_inter = {pos, end};
                interval result;
                sv_consensus cons;
                deletion(contig_index, del_begin, del_end, sv_inter, targs, &result, &cons);

                pthread_mutex_lock(targs->out_err_mutex);
                printf("(DEL) contig: %d, org pos: %u, org end: %u, ref pos: ", contig_index, sv_inter.start, sv_inter.end);
                if (result.start == 0xFFFFFFFF) printf("NA, ref end: "); else printf("%d, ref end: ", result.start);
                if (result.end == 0xFFFFFFFF) printf("NA, "); else printf("%d, ", result.end);

                if (result.start == 0xFFFFFFFF) printf("diff pos: NA, "); else printf("diff pos: %d, ", result.start - pos);
                if (result.end == 0xFFFFFFFF) printf("diff end: NA"); else printf("diff end: %d", result.end - end);
                if (cons.seq != NULL) printf(", cons_len: %ld, support: %ld\n", cons.len, cons.support); else printf(", cons_len: NA\n");
                pthread_mutex_unlock(targs->out_err_mutex);

                sv_consensus_free(&cons);
            }
            break;
        case SV_INV:
            {
                svtrek_index inv_begin_start = pos > (uint32_t)targs->wider_interval ? pos - targs->wider_interval : 1;
                svtrek_index inv_end_start = end > (uint32_t)targs->wider_interval ? end - targs->wider_interval : 1;
                interval inv_begin = {inv_begin_start, pos + targs->wider_interval};
                interval inv_end = {inv_end_start, end + targs->wider_interval};
                interval sv_inter = {pos, end};
                interval result;
                inversion(contig_index, inv_begin, inv_end, sv_inter, targs, &result);

                pthread_mutex_lock(targs->out_err_mutex);
                printf("(INV) chr: %d, org pos: %u, org end: %u, ref pos: ", contig_index, pos, end);
                if (result.start == 0xFFFFFFFF) printf("NA, ref end: "); else printf("%u, ref end: ", result.start);
                if (result.end == 0xFFFFFFFF) printf("NA\n"); else printf("%u\n", result.end);
                pthread_mutex_unlock(targs->out_err_mutex);
            }
            break;
        default:
#ifdef DEBUG
            fprintf(stderr, "[FAIL] Variant %s, begin: %u, end: %u, %s Unkown type.\n", contig, pos, end, sv_to_str(sv_type));
#endif
            break;
        }

        free(line);
    }

    sam_close(targs->hargs.fp_in);
    hts_idx_destroy(targs->hargs.bam_file_index);
    bam_hdr_destroy(targs->hargs.bam_hdr);

    targs->hargs.fp_in = NULL;
    targs->hargs.bam_hdr = NULL;
    targs->hargs.bam_file_index = NULL;
}

int process_vcf(audt_args *params) {

    printf("[INFO] Started processing variation file.\n");

    t_arg *t_args = (t_arg*)malloc(params->thread_number * sizeof(t_arg));
    pthread_mutex_t queue_mutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_t out_err_mutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_cond_t cond_not_full = PTHREAD_COND_INITIALIZER;
    pthread_cond_t cond_not_empty = PTHREAD_COND_INITIALIZER;
    int exit_signal = 0;

    pthread_mutex_init(&queue_mutex, NULL);
    pthread_mutex_init(&out_err_mutex, NULL);
    pthread_cond_init(&cond_not_full, NULL);
    pthread_cond_init(&cond_not_empty, NULL);

    line_queue queue;
    line_queue_init(&queue, params->tload_factor * params->thread_number);

    for (int i=0; i<params->thread_number; i++) {
        t_args[i].hargs.fp_in = hts_open(params->bam_file, "r");
        t_args[i].hargs.bam_hdr = sam_hdr_read(t_args[i].hargs.fp_in);
        t_args[i].hargs.bam_file_index = sam_index_load(t_args[i].hargs.fp_in, params->bam_file);
        t_args[i].wider_interval = params->wider_interval;
        t_args[i].median_interval = params->median_interval;
        t_args[i].narrow_interval = params->narrow_interval;
        t_args[i].consensus_interval_range = params->consensus_interval_range;
        t_args[i].consensus_interval = params->consensus_interval;
        t_args[i].consensus_min_count = params->consensus_min_count;
        t_args[i].queue = &(queue);
        t_args[i].queue_mutex = &queue_mutex;
        t_args[i].out_err_mutex = &out_err_mutex;
        t_args[i].cond_not_full = &cond_not_full;
        t_args[i].cond_not_empty = &cond_not_empty;
        t_args[i].exit_signal = &exit_signal;
    }

    struct tpool *tm;

    tm = tpool_create(params->thread_number);

    for (int i=0; i<params->thread_number; i++) {
        tpool_add_work(tm, thread_func, t_args+i);
    }

    uint64_t current_size = 1048576;
    char *line = (char *)malloc(current_size);

    FILE *file = fopen(params->vcf_file, "r");

    int line_count = 0;
    while (fgets(line, current_size, file) != NULL) {
        line_count++;
        size_t len = strlen(line);
        int skip_line = 0;

        while (len == current_size - 1 && line[len - 1] != '\n') {
            current_size *= 2;
            char *temp_line = (char *)realloc(line, current_size);
            if (!temp_line) {
                fprintf(stderr, "VCF: Memory reallocation failed.\n");
                free(line);
                fclose(file);
                return -1;
            }
            line = temp_line;

            if (fgets(line + len, current_size - len, file) == NULL) {
                skip_line = 1;
                break;
            }
            len = strlen(line);
        }

        if (skip_line || len < 2 || line[0] == '#')
            continue;

        if (line[len - 1] == '\n') {
            line[len - 1] = '\0';
            len--;
        }

        char *queue_line = strdup(line);
        if (!queue_line) {
            fprintf(stderr, "VCF: Memory allocation failed.\n");
            continue;
        }
        line_queue_push(&queue, queue_line, &queue_mutex, &cond_not_full, &cond_not_empty);
    }

    fclose(file);

    // Set the exit flag and wake workers while holding queue_mutex. Workers
    // read exit_signal under queue_mutex inside line_queue_pop(); setting it
    // without the lock risks a lost wakeup (a worker checks exit_signal == 0,
    // then blocks in pthread_cond_wait after this broadcast has already fired).
    pthread_mutex_lock(&queue_mutex);
    exit_signal = 1;
    pthread_cond_broadcast(&cond_not_empty);
    pthread_mutex_unlock(&queue_mutex);
    tpool_wait(tm);
    tpool_destroy(tm);
    pthread_mutex_destroy(&queue_mutex);
    pthread_mutex_destroy(&out_err_mutex);
    pthread_cond_destroy(&cond_not_full);
    pthread_cond_destroy(&cond_not_empty);

    free(queue.lines);
    free(t_args);

    printf("[INFO] Ended processing variation file\n");

    return 0;
}

int audit(int argc, char *argv[]) {

    audt_args params;
    init_audt(argc, argv, &params);

    process_vcf(&params);

    return 1;
}