#include "audit.h"

sv_type_t parse_sv_type(const char *sv_str) {
    if (strcmp(sv_str, "INS") == 0 || strcmp(sv_str, "INS:ME") == 0 ) return SV_INS;
    if (strcmp(sv_str, "DEL") == 0 || strcmp(sv_str, "DEL:ME") == 0 ) return SV_DEL;
    if (strcmp(sv_str, "INV") == 0) return SV_INV;
    if (strcmp(sv_str, "DUP") == 0) return SV_DUP;
    if (strcmp(sv_str, "TRA") == 0) return SV_TRA;
    if (strcmp(sv_str, "BND") == 0) return SV_BND;
    return SV_UNKNOWN;
}

void line_queue_init(line_queue *queue, int capacity) {
    queue->lines = (char **)malloc(capacity * sizeof(char *));
    queue->size = 0;
    queue->capacity = capacity;
    queue->front = 0;
    queue->rear = 0;
}

void line_queue_push(line_queue *queue, char *line, pthread_mutex_t *queue_mutex, pthread_cond_t *cond_not_full, pthread_cond_t *cond_not_empty) {
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

char *line_queue_pop(line_queue *queue, pthread_mutex_t *queue_mutex, pthread_cond_t *cond_not_full, pthread_cond_t *cond_not_empty, int *exit_signal) {
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

void thread_func(void *_params) {

    t_arg *targs = (t_arg *)_params;
    line_queue *queue = targs->queue;

    while (1) {
        // split the line by tab characters
        char *line = line_queue_pop(queue, targs->queue_mutex, targs->cond_not_full, targs->cond_not_empty, targs->exit_signal);
        if (line == NULL) {
            break;
        }

        char *chrom, *index, *seq, *alt, *info;
        char *saveptr;

        chrom = strtok_r(line, "\t", &saveptr); // get chromosome name
        index = strtok_r(NULL, "\t", &saveptr); // get index  
        if (index == NULL) {
            fprintf(stderr, "VCF: no index at line: %s\n", line);
            free(line);
            continue;
        }
        
        strtok_r(NULL, "\t", &saveptr); // skip ID
        seq = strtok_r(NULL, "\t", &saveptr); // skip sequence
        int seq_len = strlen(seq);

        alt = strtok_r(NULL, "\t", &saveptr); // skip ALT alleles
        char *alt_saveptr;
        char *alt_token = strtok_r(alt, ",", &alt_saveptr); // split ALT alleles by comma
        size_t max_alt_len = 0;
        size_t min_alt_len = 0x7FFFFFFF;
        while (alt_token != NULL) {
            if (strlen(alt_token) > max_alt_len) {
                max_alt_len = strlen(alt_token);
            }
            if (strlen(alt_token) < min_alt_len) {
                min_alt_len = strlen(alt_token);
            }
            alt_token = strtok_r(NULL, ",", &alt_saveptr);
        }

        strtok_r(NULL, "\t", &saveptr); // skip QUAL alleles
        strtok_r(NULL, "\t", &saveptr); // skip FILTER alleles
        info = strtok_r(NULL, "\t", &saveptr); // get INFO alleles
        
        int chrom_index;
        uint32_t pos, end;
        sv_type_t sv_type_enum = SV_UNKNOWN;

        // 1. Parse CHROM
        if (strncmp(chrom, "chr", 3) == 0) {
            chrom_index = atoi(chrom + 3);  // Skip "chr"
        } else {
            chrom_index = atoi(chrom);
        }

        // 2. Parse POS
        pos = strtol(index, NULL, 10);
        if (pos == 0 && index[0] != '0') { // Check for conversion error
            fprintf(stderr, "[ERROR] Conversion error to pos %s\n", index);
            free(line);
            continue;
        }

        // 3. Parse SVTYPE=
        char *sv_start = strstr(info, "SVTYPE=");
        if (sv_start) {
            sv_start += 7; // move past "SVTYPE="
            char *sv_end = strchr(sv_start, ';');
            if (!sv_end) {
                sv_end = sv_start + strlen(sv_start);
            }
            char sv_buf[16] = {0};
            size_t len = sv_end - sv_start;
            if (len >= sizeof(sv_buf)) {
                len = sizeof(sv_buf) - 1;
            }
            strncpy(sv_buf, sv_start, len);
            sv_buf[len] = '\0';

            sv_type_enum = parse_sv_type(sv_buf);
        } else {
            if (seq_len == 1 && __SV_MIN_LENGTH < max_alt_len) {
                sv_type_enum = SV_INS;
            } else if (__SV_MIN_LENGTH < seq_len && min_alt_len == 1) {
                sv_type_enum = SV_DEL;
            } else {
                sv_type_enum = SV_UNKNOWN;
                free(line);
                continue;
            }
        }

        // 4. Parse END=
        char *end_start = strstr(info, "END=");
        if (end_start) {
            end_start += 4; // move past "END="
            char *end_end = strchr(end_start, ';');
            if (!end_end) {
                end_end = end_start + strlen(end_start);
            }
            char end_buf[32] = {0};
            size_t len = end_end - end_start;
            if (len >= sizeof(end_buf)) {
                len = sizeof(end_buf) - 1;
            }
            strncpy(end_buf, end_start, len);
            end_buf[len] = '\0';
            end = strtol(end_buf, NULL, 10);
            if (end == 0 && end_buf[0] != '0') {
                free(line);
                continue;
            }
        } else {
            end = pos + strlen(seq);
        }

        if (sv_type_enum == SV_DEL || sv_type_enum == SV_INV) {
            if (end - pos < __SV_MIN_LENGTH) {
                free(line);
                continue;
            }
        }

        switch (sv_type_enum) {
        case SV_INS:
            {
                uint32_t begin_start = pos > (uint32_t)targs->median_interval ? pos - targs->median_interval : 1;
                interval begin = {begin_start, pos + targs->median_interval};
                uint32_t result;
                sv_consensus cons;
                insertion(chrom_index, begin, pos, targs, &result, &cons);

                pthread_mutex_lock(targs->out_err_mutex);
                if (result == 0xFFFFFFFF) {
                    printf("(INS) chr: %d, org pos: %u, ref pos: NA", chrom_index, pos);
                } else {
                    printf("(INS) chr: %d, org pos: %u, ref pos: %u, diff: %d", chrom_index, pos, result, result-pos);
                }
                if (cons.seq != NULL) {
                    printf(", cons_len: %d, support: %d, cons: %.60s%s\n", cons.len, cons.support, cons.seq, cons.len > 60 ? "..." : "");
                } else {
                    printf(", cons_len: NA\n");
                }
                pthread_mutex_unlock(targs->out_err_mutex);

                sv_consensus_free(&cons);
            }
            break;
        case SV_DEL:
            {
                if (__SV_MIN_LENGTH < end-pos) {
                    uint32_t del_begin_start = pos > (uint32_t)targs->wider_interval ? pos - targs->wider_interval : 1;
                    uint32_t del_end_start = end > (uint32_t)targs->narrow_interval ? end - targs->narrow_interval : 1;
                    interval del_begin = {del_begin_start, pos + targs->narrow_interval};
                    interval del_end = {del_end_start, end + targs->narrow_interval};
                    interval sv_inter = {pos, end};
                    interval result;
                    sv_consensus cons;
                    deletion(chrom_index, del_begin, del_end, sv_inter, targs, &result, &cons);

                    pthread_mutex_lock(targs->out_err_mutex);
                    printf("(DEL) chr: %d, org pos: %u, org end: %u, ref pos: ", chrom_index, sv_inter.start, sv_inter.end);
                    if (result.start == 0xFFFFFFFF) {
                        printf("NA, ref end: ");
                    } else {
                        printf("%d, ref end: ", result.start);
                    }
                    if (result.end == 0xFFFFFFFF) {
                        printf("NA, ");
                    } else {
                        printf("%d, ", result.end);
                    }

                    if (result.start == 0xFFFFFFFF) {
                        printf("diff pos: NA, ");
                    } else {
                        printf("diff pos: %d, ", result.start-pos);
                    }
                    if (result.end == 0xFFFFFFFF) {
                        printf("diff end: NA");
                    } else {
                        printf("diff end: %d", result.end-end);
                    }
                    if (cons.seq != NULL) {
                        printf(", cons_len: %d, support: %d\n", cons.len, cons.support);
                    } else {
                        printf(", cons_len: NA\n");
                    }
                    pthread_mutex_unlock(targs->out_err_mutex);

                    sv_consensus_free(&cons);
                }
            }
            break;
        case SV_INV:
            {
                if (__SV_MIN_LENGTH < end-pos) {
                    uint32_t inv_begin_start = pos > (uint32_t)targs->wider_interval ? pos - targs->wider_interval : 1;
                    uint32_t inv_end_start = end > (uint32_t)targs->wider_interval ? end - targs->wider_interval : 1;
                    interval inv_begin = {inv_begin_start, pos + targs->wider_interval};
                    interval inv_end = {inv_end_start, end + targs->wider_interval};
                    interval sv_inter = {pos, end};
                    interval result;
                    inversion(chrom_index, inv_begin, inv_end, sv_inter, targs, &result);

                    pthread_mutex_lock(targs->out_err_mutex);
                    printf("(INV) chr: %d, org pos: %u, org end: %u, ref pos: ", chrom_index, pos, end);
                    if (result.start == 0xFFFFFFFF) printf("NA, ref end: "); else printf("%u, ref end: ", result.start);
                    if (result.end == 0xFFFFFFFF) printf("NA\n"); else printf("%u\n", result.end);
                    pthread_mutex_unlock(targs->out_err_mutex);
                }
            }
            break;
        default:
            fprintf(stderr, "[ERROR] Unkown type.\n");
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

    exit_signal = 1;
    pthread_cond_broadcast(&cond_not_empty);
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