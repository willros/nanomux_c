#define NOB_IMPLEMENTATION
#include "nob.h"
#include <zlib.h>
#include "kseq.h"
#include <limits.h> 


uint8_t* edit_distance(const uint8_t *string, size_t stringlen, const uint8_t *pat, size_t patlen, int k) {
    size_t len1 = stringlen + 1;
    size_t len2 = patlen + 1;

    unsigned int* dist_holder = (unsigned int*)malloc(len1 * len2 * sizeof(unsigned int));
    if (!dist_holder) {
        return NULL; 
        nob_log(NOB_ERROR, "Buy more ram lol");
        exit(1);
    }

    dist_holder[0] = 0;
    for (size_t j = 1; j < len2; ++j) {
        dist_holder[j] = j;
    }
    for (size_t i = 1; i < len1; ++i) {
        dist_holder[i * len2] = 0;
    }

    int best = INT_MAX;
    unsigned int end = 0;

    for (size_t j = 1; j < len2; ++j) {
        int any_below_threshold = 0; // Flag used for early exit
        for (size_t i = 1; i < len1; ++i) {
            int sub = (string[i - 1] == pat[j - 1]) ? 0 : 1;

            if (sub == 0) {
                dist_holder[i * len2 + j] = dist_holder[(i - 1) * len2 + (j - 1)];
            } else {
                dist_holder[i * len2 + j] = dist_holder[(i - 1) * len2 + j] + 1;
                if (dist_holder[i * len2 + j] > dist_holder[i * len2 + (j - 1)] + 1) {
                    dist_holder[i * len2 + j] = dist_holder[i * len2 + (j - 1)] + 1;
                }
                if (dist_holder[i * len2 + j] > dist_holder[(i - 1) * len2 + (j - 1)] + 1) {
                    dist_holder[i * len2 + j] = dist_holder[(i - 1) * len2 + (j - 1)] + 1;
                }
            }

            if (dist_holder[i * len2 + j] <= k) {
                any_below_threshold = 1;
            }

            if (j == (len2 - 1) && dist_holder[i * len2 + j] <= k && dist_holder[i * len2 + j] < best) {
                best = dist_holder[i * len2 + j];
                end = i;
            }
        }

        if (!any_below_threshold) {
            free(dist_holder);
            return NULL; 
        }
    }

    free(dist_holder);

    if (best <= k) {
        return (uint8_t*)(string + end - 1); // Return the end position of alignment
    } else {
        return NULL;
    }
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
    size_t length;
} Read;

typedef struct {
    Read *items;
    size_t count;
    size_t capacity;
} Reads;


KSEQ_INIT(gzFile, gzread)
Reads parse_fastq(const char *fastq_file_path, int read_len_min, int read_len_max) {
    
    gzFile fp = gzopen(fastq_file_path, "r"); 
    kseq_t *seq = kseq_init(fp); 
    int l;
    Reads reads = {0};
    
    while ((l = kseq_read(seq)) >= 0) { 
        size_t length = strlen(seq->seq.s);
        
        if (length >= read_len_min && length <= read_len_max) {
            Read read = {0};
            read.seq = (const uint8_t *)strdup(seq->seq.s);
            read.name = strdup(seq->name.s);
            read.qual = strdup(seq->qual.s);
            read.length = length;
            nob_da_append(&reads, read);
        }
    }
    
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
        fprintf(stderr, "Error: Memory allocation failed for reverse complement sequence.\n");
        exit(1); 
    }

    complement_sequence(barcode->rv, comp_seq, barcode->rv_length);
    barcode->rv_comp = comp_seq;
}

void compute_reverse_complement_fw(Barcode *barcode) {
    uint8_t *comp_seq = malloc((barcode->fw_length + 1) * sizeof(uint8_t));
    if (comp_seq == NULL) {
        fprintf(stderr, "Error: Memory allocation failed for reverse complement sequence.\n");
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


void print_barcode(Barcode *b) {
    printf("FW: %.*s, RV: %.*s\n",
           (int)b->fw_length, b->fw,
           (int)b->rv_length, b->rv);
}

void print_reads(Read *r) {
    printf("Read: %.*s\n", (int)r->length, r->seq);
}


void reverse_complement_barcode_rv(Barcode *barcode, uint8_t *comp_seq) {
    complement_sequence(barcode->rv, comp_seq, barcode->rv_length);
}


void substring(const uint8_t *source, size_t start, size_t length, uint8_t *dest) {
    memcpy(dest, source + start, length);
    dest[length] = '\0';
}


void append_read_to_gzip_fastq(gzFile gzfp, Read *read) {
    gzprintf(gzfp, "@%s\n", read->name);
    gzprintf(gzfp, "%s\n", read->seq);
    gzprintf(gzfp, "+\n");
    gzprintf(gzfp, "%s\n", read->qual);
}


int main(int argc, char **argv) {    
    
    const char *program = nob_shift_args(&argc, &argv);
    if (argc < 7) {
        nob_log(NOB_ERROR, "usage: %s <barcode_file: path> <fastq_file: path> <read_len_min: int> <read_len_max: int> <barcode_pos: int> <k: int> <output folder: path>", program);
        return 1;
    }
    
    const char *barcode_file = nob_shift_args(&argc, &argv);
    const char *fastq_file = nob_shift_args(&argc, &argv);
    int read_len_min = atoi(nob_shift_args(&argc, &argv));
    int read_len_max = atoi(nob_shift_args(&argc, &argv));
    size_t barcode_pos = atoi(nob_shift_args(&argc, &argv));
    int k = atoi(nob_shift_args(&argc, &argv));
    const char *output = nob_shift_args(&argc, &argv);
    nob_log(NOB_INFO, "Read len min: %i", read_len_min);
    nob_log(NOB_INFO, "Read len max: %i", read_len_max);
    nob_log(NOB_INFO, "Barcode position: 0 - %zu", barcode_pos);
    nob_log(NOB_INFO, "k: %i", k);
    nob_log(NOB_INFO, "Output folder: %s", output);

    if (!nob_mkdir_if_not_exists(output)) return 1;
    
    
    // create log file
    char log_file[256];
    snprintf(log_file, sizeof(log_file), "%s/nanomux.log", output);
    FILE *LOG_FILE = fopen(log_file, "ab");
    
    fprintf(LOG_FILE, "Nanomux\n\n");
    fprintf(LOG_FILE, "Barcodes: %s\n", barcode_file);
    fprintf(LOG_FILE, "Fastq: %s\n", fastq_file);
    fprintf(LOG_FILE, "Read len min: %i\n", read_len_min);
    fprintf(LOG_FILE, "Read len max: %i\n", read_len_max);
    fprintf(LOG_FILE, "Barcode position: %zu\n", barcode_pos);
    fprintf(LOG_FILE, "K: %i\n", k);
    fprintf(LOG_FILE, "Output folder: %s\n", output);
    
    nob_log(NOB_INFO, "Parsing barcode file %s", barcode_file);
    Nob_String_Builder sb = {0};
    if (!nob_read_entire_file(barcode_file, &sb)) return 1;
    Nob_String_View content = nob_sv_from_parts(sb.items, sb.count);
    Barcodes barcodes = parse_barcodes(content);
    
    nob_log(NOB_INFO, "Parsing fastq file %s", fastq_file);
    Reads reads = parse_fastq(fastq_file, read_len_min, read_len_max);
    nob_log(NOB_INFO, "Number of reads: %zu", reads.count);
    fprintf(LOG_FILE, "Number of reads: %zu\n", reads.count);
    
    fclose(LOG_FILE);
    
    uint8_t first_part[1000];
    uint8_t last_part[1000];
    
    for (size_t i = 0; i < barcodes.count; ++i) {
        
        int counter = 0;
        Barcode b = barcodes.items[i];
        nob_log(NOB_INFO, "Name: "SV_Fmt, SV_Arg(b.name));
        nob_log(NOB_INFO, "Fw: %.*s", b.fw_length, b.fw);
        nob_log(NOB_INFO, "rv comp: %.*s", b.rv_length, b.rv_comp);
        
        // save fastq
        char fastq_name[256]; 
        snprintf(fastq_name, sizeof(fastq_name), "%s/%.*s.fq.gz", output, (int)b.name.count, b.name.data);
        
        gzFile new_fastq = gzopen(fastq_name, "ab");
        if (!new_fastq) {
            nob_log(NOB_INFO, "Failed to open new fastq_file for appending");
            return 1;
        }
        
        for (size_t j = 0; j < reads.count; ++j) {
            
            Read r = reads.items[j];
            substring(r.seq, 0, barcode_pos, first_part); 
            size_t last_part_start = r.length - barcode_pos;
            substring(r.seq, last_part_start, barcode_pos, last_part); 
            
            uint8_t *match_first_fw = edit_distance(first_part, barcode_pos, b.fw, b.fw_length, k);
            
            if (match_first_fw != NULL) {
                uint8_t *match_last_fw = edit_distance(last_part, barcode_pos, b.rv_comp, b.rv_length, k);
                if (match_last_fw != NULL) {
                    counter ++;
                    append_read_to_gzip_fastq(new_fastq, &r);
                    continue;
                }
            }
            
            uint8_t *match_first_rv = edit_distance(first_part, barcode_pos, b.rv, b.rv_length, k);
            if (match_first_rv != NULL) {
                uint8_t *match_last_rv = edit_distance(last_part, barcode_pos, b.fw_comp, b.fw_length, k);
                if (match_last_rv != NULL) {
                    counter ++;
                    append_read_to_gzip_fastq(new_fastq, &r);
                }
            }
        }
        nob_log(NOB_INFO, "matches: %i", counter);
        nob_log(NOB_INFO, "Saving fastq file: %s", fastq_name);
        printf("\n");
        gzclose(new_fastq);
    }
    
    nob_log(NOB_INFO, "Done!");
    nob_da_free(barcodes);
    nob_da_free(reads);
    return 0;
}
