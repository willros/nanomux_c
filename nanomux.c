//TODO: readsplitting
// this seems to be the adapter: TACTTCGTTCAGTTACGTATTGCT
// they also use just the first 14 bp: TACTTCGTTCAGTT
// search for that, strstr()? and just split by that into two reads? 

// Copyright 2024 William Rosenbaum <william.rosenbaum88@gmail.com>
//
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#define NOB_IMPLEMENTATION
#include "nob.h"
#include <zlib.h>
#include "kseq.h"
#include <limits.h> 
#include <stdint.h>
#include <pthread.h>


// inspired by flexplex: https://github.com/DavidsonGroup/flexiplex
int edit_distance(const uint8_t* text, size_t textlen, const uint8_t* pattern, size_t patlen, int k) {
    int len1 = textlen + 1;
    int len2 = patlen + 1;

    int* dist_holder = (int*)malloc(len1 * len2 * sizeof(int));
    if (!dist_holder) return -1; 

    dist_holder[0] = 0; //[0][0]
    for (int j = 1; j < len2; ++j) dist_holder[j] = j; //[0][j]
    for (int i = 1; i < len1; ++i) dist_holder[i * len2] = 0; //[i][0]

    int best = INT_MAX; 
    int end = len1 - 1;

    for (int j = 1; j < len2; ++j) {
        int any_below_threshold = 0; 
        for (int i = 1; i < len1; ++i) {
            int sub = (text[i - 1] == pattern[j - 1]) ? 0 : 1; 

            if (sub == 0) {
                dist_holder[i * len2 + j] = dist_holder[(i - 1) * len2 + (j - 1)];
            } else {
                dist_holder[i * len2 + j] = dist_holder[(i - 1) * len2 + j] + 1; // [i-1][j]
                if (dist_holder[i * len2 + j] > dist_holder[i * len2 + (j - 1)] + 1) // [i][j-1]
                    dist_holder[i * len2 + j] = dist_holder[i * len2 + (j - 1)] + 1;
                if (dist_holder[i * len2 + j] > dist_holder[(i - 1) * len2 + (j - 1)] + 1) // [i-1][j-1]
                    dist_holder[i * len2 + j] = dist_holder[(i - 1) * len2 + (j - 1)] + 1;
            }
            if (dist_holder[i * len2 + j] <= k) {
                any_below_threshold = 1;
            }

            if (j == (len2 - 1) && dist_holder[i * len2 + j] < best) {
                best = dist_holder[i * len2 + j];
                end = i; // Update the end position of alignment
                if (best <= k) {
                    free(dist_holder);
                    return end; 
                }
            }
        }

        if (!any_below_threshold) {
            free(dist_holder);
            return -1; 
        }
    }

    free(dist_holder);
    return -1; 
}

typedef enum {
    BARCODE_NAME = 0,
    BARCODE_FW,
    BARCODE_RV
} Barcode_Attr;

typedef struct {
    Nob_String_View name;
    const uint8_t *fw;
    size_t fw_length;
    const uint8_t *rv;
    size_t rv_length;
    const uint8_t *rv_comp;  
    const uint8_t *fw_comp;  
} Barcode;

typedef struct {
    Barcode *items;
    size_t count;
    size_t capacity;
} Barcodes;

typedef struct {
    const uint8_t *seq;
    const char *name;
    const char *qual;
    int length;
} Read;

void free_read(Read *read) {
    if (read) {
        free((void*)read->seq); 
        free((void*)read->qual); 
        free((void*)read->name); 
    }
}

typedef struct {
    Read *items;
    size_t count;
    size_t capacity;
} Reads;


typedef struct {
    Barcode barcode;
    Reads reads;
    const char *output;
    size_t barcode_pos;
    int k;
    int trim;
    FILE *S_FILE;
    pthread_mutex_t *s_mutex;
} ThreadData;


// TODO: update the adapter len when changing the adapter
// TODO: change the parameters of the adapter position and so on
const uint8_t *ADAPTER = (const uint8_t*)"TACTTCGTTCAGTTACGTATTGCT";
#define ADAPTER_LEN 24
#define CONCATENATED_READ_LEN_MIN 2500
#define CONCATENATED_READ_LEN_MAX 5000
#define ADAPTER_POSITION 1200

void split_reads_by_adapter(const Read *read, Reads *reads, int *counter, int read_len_min, int read_len_max) {

    int match = edit_distance(read->seq, read->length, ADAPTER, ADAPTER_LEN, 1); // k = 1
    
    if (match > ADAPTER_POSITION) {  
        *counter += 1;
        // two loops, one for each new read
        for (int i = 1; i < 3; i++) {
            
            int start = (i == 1) ? 0 : match;
            int end = (i == 1) ? match : read->length;
            int new_length = end - start;
            
            if (new_length >= read_len_min && new_length <= read_len_max) {
                char new_read_name[100];
                sprintf(new_read_name, "%s_%i", read->name, i);

                char new_read_seq[3500];
                snprintf(new_read_seq, sizeof(new_read_seq), "%.*s", (int)new_length, read->seq + start);

                char new_read_qual[3500];
                snprintf(new_read_qual, sizeof(new_read_qual), "%.*s", (int)new_length, read->qual + start);

                // add everything to the da Reads
                Read new_read = {0};
                new_read.seq = (const uint8_t *)strdup(new_read_seq);
                new_read.name = strdup(new_read_name);
                new_read.qual = strdup(new_read_qual);
                new_read.length = new_length;
                nob_da_append(reads, new_read);
            }
        }
    }
}

KSEQ_INIT(gzFile, gzread)
Reads parse_fastq(const char *fastq_file_path, int read_len_min, int read_len_max, bool split_reads) {
    
    gzFile fp = gzopen(fastq_file_path, "r"); 
    if (!fp) {
        nob_log(NOB_INFO, "Failed to open fastq file, exiting");
        exit(1);
    }
    
    kseq_t *seq = kseq_init(fp); 
    int l;
    Reads reads = {0};
    
    int counter = 0;
    
    while ((l = kseq_read(seq)) >= 0) { 
        int length = strlen(seq->seq.s);
        
        if (split_reads) {
            if (length >= CONCATENATED_READ_LEN_MIN && length <= CONCATENATED_READ_LEN_MAX) {
                Read concated_read = {0};
                concated_read.seq = (const uint8_t *)strdup(seq->seq.s);
                concated_read.name = strdup(seq->name.s);
                concated_read.qual = strdup(seq->qual.s);
                concated_read.length = length;

                split_reads_by_adapter(&concated_read, &reads, &counter, read_len_min, read_len_max);
                free_read(&concated_read);
                continue;
            }
        }
        
        if (length >= read_len_min && length <= read_len_max) {
            Read read = {0};
            read.seq = (const uint8_t *)strdup(seq->seq.s);
            read.name = strdup(seq->name.s);
            read.qual = strdup(seq->qual.s);
            read.length = length;
            nob_da_append(&reads, read);
        }
    }
    
    if (split_reads) nob_log(NOB_INFO, "Number of splitted reads: %i", counter);
    
    kseq_destroy(seq); 
    gzclose(fp); 
    return reads;
}

uint8_t complement(uint8_t nucleotide) {
    switch (nucleotide) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        default: return 'N'; 
    }
}

void complement_sequence(const uint8_t *seq, uint8_t *comp_seq, size_t length) {
    for (size_t i = 0; i < length; i++) {
        comp_seq[length - 1 - i] = complement(seq[i]);
    }
    comp_seq[length] = '\0';
}

void compute_reverse_complement_rv(Barcode *barcode) {
    uint8_t *comp_seq = malloc((barcode->rv_length + 1) * sizeof(uint8_t));
    if (comp_seq == NULL) {
        nob_log(NOB_ERROR, "Memory allocation failed for reverse complement, exiting");
        exit(1); 
    }

    complement_sequence(barcode->rv, comp_seq, barcode->rv_length);
    barcode->rv_comp = comp_seq;
}

void compute_reverse_complement_fw(Barcode *barcode) {
    uint8_t *comp_seq = malloc((barcode->fw_length + 1) * sizeof(uint8_t));
    if (comp_seq == NULL) {
        nob_log(NOB_ERROR, "Memory allocation failed for reverse complement, exiting");
        exit(1); 
    }

    complement_sequence(barcode->fw, comp_seq, barcode->fw_length);
    barcode->fw_comp = comp_seq;
}

Barcodes parse_barcodes(Nob_String_View content) {
    Barcodes barcodes = {0};

    while (content.count > 0) {
        Nob_String_View line = nob_sv_chop_by_delim(&content, '\n');
        Barcode barcode = {0};

        for (int i = 0; line.count > 0; ++i) {
            Nob_String_View attr = nob_sv_chop_by_delim(&line, ',');

            switch (i) {
                case BARCODE_NAME:
                    barcode.name = attr;
                    break;
                case BARCODE_FW:
                    barcode.fw = (const uint8_t *)attr.data;
                    barcode.fw_length = attr.count;
                    compute_reverse_complement_fw(&barcode); 
                    break;
                case BARCODE_RV:
                    barcode.rv = (const uint8_t *)attr.data;
                    barcode.rv_length = attr.count;
                    compute_reverse_complement_rv(&barcode); 
                    break;
                default:
                    break;
            }
        }
        nob_da_append(&barcodes, barcode);
    }

    return barcodes;
}

void substring(const uint8_t *source, size_t start, size_t length, uint8_t *dest) {
    memcpy(dest, source + start, length);
    dest[length] = '\0';
}

void append_read_to_gzip_fastq_trim(gzFile gzfp, Read *read, int start, int end) {
    if (start < 0) start = 0;
    if (end > (int)read->length) end = read->length;
    size_t trimmed_length = end - start;
    gzprintf(gzfp, "@%s\n", read->name);
    gzprintf(gzfp, "%.*s\n", (int)trimmed_length, read->seq + start);
    gzprintf(gzfp, "+\n");
    gzprintf(gzfp, "%.*s\n", (int)trimmed_length, read->qual + start);
}

void append_read_to_gzip_fastq(gzFile gzfp, Read *read) {
    gzprintf(gzfp, "@%s\n", read->name);
    gzprintf(gzfp, "%s\n", read->seq);
    gzprintf(gzfp, "+\n");
    gzprintf(gzfp, "%s\n", read->qual);
}

void process_single_barcode_thread(
    const Barcode b,
    const Reads reads,
    const char *output,
    const size_t barcode_pos,
    const int k,
    const int trim,
    FILE *S_FILE,
    pthread_mutex_t *s_mutex
) {
    int counter = 0;
    
    // Save fastq
    char fastq_name[256]; 
    snprintf(fastq_name, sizeof(fastq_name), "%s/%.*s.fq.gz", output, (int)b.name.count, b.name.data);
    
    gzFile new_fastq = gzopen(fastq_name, "ab");
    if (!new_fastq) {
        nob_log(NOB_INFO, "Failed to open new fastq_file for appending, exiting");
        return;
    }
    
    // TODO: change this to read_len_max
    uint8_t first_part[1000];
    uint8_t last_part[1000];
    
    for (size_t j = 0; j < reads.count; ++j) {
        int match_first_fw = -1;
        int match_last_fw = -1;
        int match_first_rv = -1;
        int match_last_rv = -1;
        
        Read r = reads.items[j];
        substring(r.seq, 0, barcode_pos, first_part); 
        size_t last_part_start = r.length - barcode_pos;
        substring(r.seq, last_part_start, barcode_pos, last_part); 
        
        match_first_fw = edit_distance(first_part, barcode_pos, b.fw, b.fw_length, k);
        // fw ------ revcomp(rv)
        if (match_first_fw != -1) {
            match_last_fw = edit_distance(last_part, barcode_pos, b.rv_comp, b.rv_length, k);
            if (match_last_fw != -1) {
                counter++;
                int trim_right = last_part_start + match_last_fw - b.rv_length;
                // trimming or not
                if (trim) {
                    append_read_to_gzip_fastq_trim(new_fastq, &r, match_first_fw, trim_right);
                } else {
                    append_read_to_gzip_fastq(new_fastq, &r);
                }
                continue;
            }
        }
        
        // rv ------ revcomp(fw)
        match_first_rv = edit_distance(first_part, barcode_pos, b.rv, b.rv_length, k);
        if (match_first_rv != -1) {
            match_last_rv = edit_distance(last_part, barcode_pos, b.fw_comp, b.fw_length, k);
            if (match_last_rv != -1) {
                int trim_right = last_part_start + match_last_rv - b.fw_length;
                counter++;
                // trimming or not
                if (trim) {
                    append_read_to_gzip_fastq_trim(new_fastq, &r, match_first_rv, trim_right);
                } else {
                    append_read_to_gzip_fastq(new_fastq, &r);
                }
            }
        }
    }
    
    nob_log(NOB_INFO, "barcode: "SV_Fmt" matches: %i", SV_Arg(b.name), counter);
    
    gzclose(new_fastq);
    
    // Write summary file with mutex lock
    pthread_mutex_lock(s_mutex);
    fprintf(S_FILE, "%.*s,%i\n", (int)b.name.count, b.name.data, counter);
    pthread_mutex_unlock(s_mutex);
}


void *thread_function(void *arg) {
    ThreadData *data = (ThreadData *)arg;
    Reads local_reads = data->reads;
    process_single_barcode_thread(
        data->barcode, 
        local_reads, 
        data->output, 
        data->barcode_pos,
        data->k, 
        data->trim,
        data->S_FILE,
        data->s_mutex
    );
    
    return NULL;
}



int main(int argc, char **argv) {    
    
    const char *program = nob_shift_args(&argc, &argv);
    
    if (argc != 10) {
        printf("[USAGE] Add the following arguments in this exact order:\n   %s\n   <barcode_file: path>\n   <fastq_file: path>\n   <read_len_min: int>\n   <read_len_max: int>\n   <barcode_pos: int>\n   <k: int>\n   <output_folder: path>\n   <trim_option: trim|notrim>\n   <split_option: split|nosplit>\n   <num_threads: int>\n", program);
        return 1;
    }
    
    const char *barcode_file = nob_shift_args(&argc, &argv);
    const char *fastq_file = nob_shift_args(&argc, &argv);
    int read_len_min = atoi(nob_shift_args(&argc, &argv));
    int read_len_max = atoi(nob_shift_args(&argc, &argv));
    size_t barcode_pos = atoi(nob_shift_args(&argc, &argv));
    int k = atoi(nob_shift_args(&argc, &argv));
    const char *output = nob_shift_args(&argc, &argv);
    const char *trim_option = nob_shift_args(&argc, &argv);
    const char *split_option = nob_shift_args(&argc, &argv);
    int num_threads = atoi(nob_shift_args(&argc, &argv));
    
    nob_log(NOB_INFO, "Read len min: %i", read_len_min);
    nob_log(NOB_INFO, "Read len max: %i", read_len_max);
    nob_log(NOB_INFO, "Barcode position: 0 -> %zu", barcode_pos);
    nob_log(NOB_INFO, "k: %i", k);
    nob_log(NOB_INFO, "Trim option: %s", trim_option);
    nob_log(NOB_INFO, "threads: %i", num_threads);
    nob_log(NOB_INFO, "Split option: %s", split_option);
    
    // handle trim option
    bool trim;
    if (strcmp(trim_option, "trim") == 0) {
        trim = true;
    } else if (strcmp(trim_option, "notrim") == 0) {
        trim = false;
    } else {
        nob_log(NOB_ERROR, "Trim option must be `trim` or `notrim`");
        return 1;
    }
    
    // handle split option
    bool split;
    if (strcmp(split_option, "split") == 0) {
        split = true;
    } else if (strcmp(split_option, "nosplit") == 0) {
        split = false;
    } else {
        nob_log(NOB_ERROR, "Split option must be `split` or `nosplit`");
        return 1;
    }

    if (!nob_mkdir_if_not_exists(output)) {
        nob_log(NOB_ERROR, "exiting");
        return 1;
    }
   
    // create log file
    char log_file[256];
    snprintf(log_file, sizeof(log_file), "%s/nanomux.log", output);
    FILE *LOG_FILE = fopen(log_file, "ab");
    if (LOG_FILE == NULL) {
        nob_log(NOB_ERROR, "Could create log file");
        return 1;
    }
    
    // create summary file
    char summary_file[256];
    snprintf(summary_file, sizeof(summary_file), "%s/nanomux_matches.csv", output);
    FILE *S_FILE = fopen(summary_file, "ab");
    if (S_FILE == NULL) {
        nob_log(NOB_ERROR, "Could not create summary file");
        return 1;
    }
    fprintf(S_FILE, "barcode,matches\n");
    
    fprintf(LOG_FILE, "Nanomux\n\n");
    fprintf(LOG_FILE, "Barcodes: %s\n", barcode_file);
    fprintf(LOG_FILE, "Fastq: %s\n", fastq_file);
    fprintf(LOG_FILE, "Read len min: %i\n", read_len_min);
    fprintf(LOG_FILE, "Read len max: %i\n", read_len_max);
    fprintf(LOG_FILE, "Barcode position: %zu\n", barcode_pos);
    fprintf(LOG_FILE, "K: %i\n", k);
    fprintf(LOG_FILE, "Output folder: %s\n", output);
    fprintf(LOG_FILE, "Trim option: %s\n", trim_option);
    fprintf(LOG_FILE, "Split option: %s\n", split_option);
    
    nob_log(NOB_INFO, "Parsing barcode file %s", barcode_file);
    Nob_String_Builder sb = {0};
    if (!nob_read_entire_file(barcode_file, &sb)) return 1;
    Nob_String_View content = nob_sv_from_parts(sb.items, sb.count);
    Barcodes barcodes = parse_barcodes(content);
    
    nob_log(NOB_INFO, "Parsing fastq file %s", fastq_file);
    Reads reads = parse_fastq(fastq_file, read_len_min, read_len_max, split);
    nob_log(NOB_INFO, "Number of reads after filtering: %zu", reads.count);
    fprintf(LOG_FILE, "Number of reads after filtering: %zu\n", reads.count);
    printf("\n");
    
    fclose(LOG_FILE);
    
    pthread_t threads[num_threads];
    ThreadData *thread_data = malloc(sizeof(ThreadData) * barcodes.count);
    if (thread_data == NULL) {
        nob_log(NOB_ERROR, "Memory allocation failed for thread data");
        return 1;
    }

    pthread_mutex_t barcode_index_mutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_t s_mutex = PTHREAD_MUTEX_INITIALIZER;
    size_t next_barcode_index = 0;

    while (next_barcode_index < barcodes.count) {
        int i = 0;
        for (; i < num_threads && next_barcode_index < barcodes.count; i++) {
            pthread_mutex_lock(&barcode_index_mutex);
            size_t current_barcode_index = next_barcode_index++;
            pthread_mutex_unlock(&barcode_index_mutex);

            thread_data[current_barcode_index] = (ThreadData){
                .barcode = barcodes.items[current_barcode_index],
                .reads = reads,
                .output = output,
                .barcode_pos = barcode_pos,
                .k = k,
                .trim = trim,
                .S_FILE = S_FILE,
                .s_mutex = &s_mutex
            };
            pthread_create(&threads[i], NULL, thread_function, &thread_data[current_barcode_index]);
        }

        for (int j = 0; j < i; j++) {
            pthread_join(threads[j], NULL);
        }
    }

    pthread_mutex_destroy(&barcode_index_mutex);
    pthread_mutex_destroy(&s_mutex);
    nob_da_free(barcodes);
    nob_da_free(reads);
    fclose(S_FILE);
    
    nob_log(NOB_INFO, "Done!");
    
    return 0;
}