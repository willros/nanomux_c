#define NOB_IMPLEMENTATION
#include "nob.h"
#include <zlib.h>
#include "kseq.h"
#include <stdint.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "thpool.h"

typedef struct {
    int min_qual;
    int min_read_len;
    int max_read_len;
    char *in_file;
    char *out_file;
    pthread_mutex_t *print_mutex;
} File; 

typedef struct {
    File *items;
    size_t count;
    size_t capacity;
} Files;

double average_qual(const char *quals, size_t len) {
    double probability_sum = 0.0;
    for (size_t i = 0; i < len; i++) {
        int phred_score = quals[i] - 33;  
        probability_sum += pow(10.0, phred_score / -10.0);
    }
    return log10(probability_sum / len) * -10.0;
}

void append_read_to_gzip_fastq(gzFile gzfp, const char *name, const char *seq, const char *qual) {
    gzprintf(gzfp, "@%s\n", name);
    gzprintf(gzfp, "%s\n", seq);
    gzprintf(gzfp, "+\n");
    gzprintf(gzfp, "%s\n", qual);
}

KSEQ_INIT(gzFile, gzread)
bool parse_fastq(
    const char *fastq_file_path, 
    int read_len_min, 
    int read_len_max, 
    int min_qual, 
    const char *out_path,
    pthread_mutex_t *print_mutex
) {
    bool result = true;
    gzFile fastq_file = gzopen(fastq_file_path, "r"); 
    if (!fastq_file) {
        nob_log(NOB_ERROR, "Failed to open %s file, exiting", fastq_file_path);
        nob_return_defer(false);
    }

    gzFile out_file = gzopen(out_path, "wb"); 
    if (!out_file) {
        nob_log(NOB_ERROR, "Failed to open %s file, exiting", out_path);
        nob_return_defer(false);
    }
    
    kseq_t *seq = kseq_init(fastq_file); 
    if (seq == NULL) {
        nob_log(NOB_ERROR, "Could not initialize %s file, exiting", fastq_file_path);
        nob_return_defer(false);
    }

    int raw_reads = 0;
    int qualified_reads = 0;
    while (kseq_read(seq) >= 0) { 
        raw_reads++;
        int length = strlen(seq->seq.s);
        if (length >= read_len_min && length <= read_len_max) {
            double average = average_qual(seq->qual.s, length);
            if (average >= min_qual) {
                qualified_reads++;
                append_read_to_gzip_fastq(
                    out_file, seq->name.s, seq->seq.s, seq->qual.s
                );
            }
        }
    }
    
    pthread_mutex_lock(print_mutex);
    nob_log(NOB_INFO, "%50s:     %15i raw reads --> %15i saved reads", fastq_file_path, raw_reads, qualified_reads);
    pthread_mutex_unlock(print_mutex);

defer:
    if (seq) kseq_destroy(seq); 
    if (fastq_file) gzclose(fastq_file); 
    if (out_file) gzclose(out_file); 
    return result;
}

void task_parse_fastq_file(void *arg) {
    File *file = (File *)arg;

    if (!parse_fastq(file->in_file, file->min_read_len, file->max_read_len, file->min_qual, file->out_file, file->print_mutex)) {
        return;
    }
}

bool is_fastq(const char *file) {
    return strstr(file, "fastq") || strstr(file, "fq");
}

bool must_be_digit(const char *arg) {
    for(; *arg != '\0'; arg++) {
        if (!isdigit(*arg)) return false;
    }
    return true;
}

char *usage = 
"[USAGE]: nanotrim -i <input> [options]\n"
"   -i    <input>             Path of folder or file\n"
"   -r    <read_length_min>   Minium length of read.    Default: 1\n"
"   -R    <read_length_max>   Minium length of read.    Default: INT_MAX\n"
"   -q    <quality>           Minimum quality of read.  Default: 1\n"
"   -t    <threads>           Number of threads to use. Default: 1\n";


int main(int argc, char **argv) {
    
    // required args
    bool i_arg = false;
    char *input;

    // optional args
    int r_arg = 1, q_arg = 1, t_arg = 1, R_arg = INT32_MAX;
    char c;
    while ((c = getopt (argc, argv, "i:r:R:q:t:")) != -1) {
        switch (c) {
            case 'i':
                i_arg = true;
                input = optarg;
                break;
            case 'r':
                if (!must_be_digit(optarg)) {
                    nob_log(NOB_ERROR, "-r must be digit");
                    nob_log(NOB_ERROR, "%s", usage);
                    return 1;
                }
                r_arg = atoi(optarg);
                break;
            case 'R':
                if (!must_be_digit(optarg)) {
                    nob_log(NOB_ERROR, "-R must be digit");
                    nob_log(NOB_ERROR, "%s", usage);
                    return 1;
                }
                R_arg = atoi(optarg);
                break;
            case 'q':
                if (!must_be_digit(optarg)) {
                    nob_log(NOB_ERROR, "-q must be digit");
                    nob_log(NOB_ERROR, "%s", usage);
                    return 1;
                }
                q_arg = atoi(optarg);
                break;
            case 't':
                if (!must_be_digit(optarg)) {
                    nob_log(NOB_ERROR, "-t must be digit");
                    nob_log(NOB_ERROR, "%s", usage);
                    return 1;
                }
                t_arg = atoi(optarg);
                break;
        }
    }
    if (!i_arg) {
        nob_log(NOB_ERROR, "%s", usage);
        return 1;
    }
    
    nob_log(NOB_INFO, "Input:               %20s", input);
    nob_log(NOB_INFO, "Minimum read length: %20i", r_arg);
    nob_log(NOB_INFO, "Maximum read length: %20i", R_arg);
    nob_log(NOB_INFO, "Minimum quality:     %20i", q_arg);
    nob_log(NOB_INFO, "Number of threads:   %20i", t_arg);

    Nob_File_Type type = nob_get_file_type(input);
    Nob_File_Paths files = {0};
    Files fastq_files = {0};
    
    pthread_mutex_t print_mutex = PTHREAD_MUTEX_INITIALIZER;

    switch (type) {
        case NOB_FILE_DIRECTORY: {
            nob_log(NOB_INFO, "`%s` is a directory", input);
            if (!nob_read_entire_dir(input, &files)) {
                nob_log(NOB_ERROR, "Failed to read directory `%s`", input);
                return 1;
            }

            char realpath[512];
            char outfile[512];
            for (size_t i = 0; i < files.count; ++i) {
                const char *file = files.items[i];

                if (strcmp(file, ".") == 0) continue;
                if (strcmp(file, "..") == 0) continue;
                if (*file == '.') continue;
                if (!is_fastq(file)) continue;

                // realpath and outfile
                snprintf(realpath, sizeof(realpath), "%s/%s", input, file);
                snprintf(outfile, sizeof(outfile), "%s/%s.filtered", input, file);

                File fastq_file = {
                    .min_qual = q_arg,
                    .min_read_len = r_arg,
                    .max_read_len = R_arg,
                    .in_file = strdup(realpath),
                    .out_file = strdup(outfile),
                    .print_mutex = &print_mutex
                };
                nob_da_append(&fastq_files, fastq_file);
            }
            break;
        }
        case NOB_FILE_REGULAR: {
            if (!is_fastq(input)) {
                nob_log(NOB_ERROR, "`%s` is not a fastq file...", input);
                return 1;
            }
            nob_log(NOB_INFO, "`%s` is a file", input);

            char outfile[512];
            snprintf(outfile, sizeof(outfile), "%s.filtered", input);

            File fastq_file = {
                .min_qual = q_arg,
                .min_read_len = r_arg,
                .max_read_len = R_arg,
                .in_file = strdup(input),
                .out_file = strdup(outfile),
                .print_mutex = &print_mutex
            };
            nob_da_append(&fastq_files, fastq_file);
            break;
        }
        default:
            nob_log(NOB_ERROR, "input: `%s` has an unknown type", input);
			return 1;
    }


    nob_log(NOB_INFO, "Generating threadpool with %i threads", t_arg);
    threadpool thpool = thpool_init(t_arg);

    for (int i = 0; i < fastq_files.count; i++) {
		thpool_add_work(thpool, task_parse_fastq_file, (void *)&fastq_files.items[i]);
    }

	thpool_wait(thpool);
	thpool_destroy(thpool);
    pthread_mutex_destroy(&print_mutex);


    nob_da_free(fastq_files);
    nob_da_free(files);

    return 0;
}
