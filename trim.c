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

#define N_THREADS 4

typedef struct {
    int min_qual;
    int min_read_len;
    char *in_file;
    char *out_file;
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
bool parse_fastq(const char *fastq_file_path, int read_len_min, int min_qual, const char *out_path) {
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
        double average = average_qual(seq->qual.s, length);
        if (length >= read_len_min && average > min_qual) {
            qualified_reads++;
            append_read_to_gzip_fastq(
                out_file, seq->name.s, seq->seq.s, seq->qual.s
            );
        }
    }
    nob_log(NOB_INFO, "%50s: %10i raw reads       <-> %10i saved reads", fastq_file_path, raw_reads, qualified_reads);

defer:
    if (seq) kseq_destroy(seq); 
    if (fastq_file) gzclose(fastq_file); 
    if (out_file) gzclose(out_file); 
    return result;
}

void task_parse_fastq_file(void *arg) {
    File *file = (File *)arg;
    printf(
        "From task thread #%lu: parsing file: %s\n", 
        (unsigned long)pthread_self(), 
        file->in_file
    );

    if (!parse_fastq(file->in_file, file->min_read_len, file->min_qual, file->out_file)) {
        return;
    }
}

bool is_fastq(const char *file) {
    return strstr(file, "fastq") || strstr(file, "fq");
}


int main(int argc, char **argv) {
    const char *program = nob_shift_args(&argc, &argv);

    if (argc != 3) {
        printf("[USAGE] Add the following arguments in this exact order:\n   %s <input> <read_len_min> <min_qual>\n", program);
        return 1;
    }

    const char *input = nob_shift_args(&argc, &argv);
    int read_len_min = atoi(nob_shift_args(&argc, &argv));
    int min_qual = atoi(nob_shift_args(&argc, &argv));


    nob_log(NOB_INFO, "input: %s", input);
    nob_log(NOB_INFO, "Read len min: %i", read_len_min);
    nob_log(NOB_INFO, "Minimum quality: %i", min_qual);

    Nob_File_Type type = nob_get_file_type(input);
    Nob_File_Paths files = {0};
    Files fastq_files = {0};

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
                    .min_qual = min_qual,
                    .min_read_len = read_len_min,
                    .in_file = strdup(realpath),
                    .out_file = strdup(outfile),
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

            char outfile[512];
            snprintf(outfile, sizeof(outfile), "%s.filtered", input);

            File fastq_file = {
                .min_qual = min_qual,
                .min_read_len = read_len_min,
                .in_file = strdup(input),
                .out_file = strdup(outfile),
            };
            nob_da_append(&fastq_files, fastq_file);
            break;
        }
        case NOB_FILE_ERROR: return 1;
        case NOB_FILE_SYMLINK: return 1;
        case NOB_FILE_OTHER: return 1;
        default:
            nob_log(NOB_ERROR, "input: `%s` has an unknown type", input);
            return 1;
    }


    nob_log(NOB_INFO, "Generating threadpool with %i threads", N_THREADS);
    threadpool thpool = thpool_init(N_THREADS);

    for (int i = 0; i < fastq_files.count; i++) {
		thpool_add_work(thpool, task_parse_fastq_file, (void *)&fastq_files.items[i]);
    }

	thpool_wait(thpool);
	thpool_destroy(thpool);

    nob_da_free(fastq_files);
    nob_da_free(files);

    return 0;
}
