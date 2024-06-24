# nanomux_c

Demultiplex your nanopore (or other) reads! Implements fuzzy matching useful for noisy reads. Written entirely in C, so should compile on most platforms.


## Quick Start
```bash
$ git clone https://github.com/willros/nanomux_c.git
$ cd nanomux_c/
$ sh build.sh
$ ./nanomux tests/bc_test.csv tests/test.fastq 600 2000 200 1 TEST_NANOMUX trim 4
```
## Help
Just run `./nanomux` to get the help message:
```bash
[ERROR] usage:
   ./nanomux <barcode_file: path> <fastq_file: path>
   <read_len_min: int> <read_len_max: int>
   <barcode_pos: int> <k: int>
   <output_folder: path>
   <trim_option: trim|notrim>
   <num_threads: int>
```
* barcode_pos: where to search for the barcode in the ends. if barcode_pos == 200 and the read length is 1000, the barcodes will be searched for from position 0 -> 200 and 800 -> 1000.
* k: allowed number of mismatches

## Usage
`barcode_file.csv` **MUST** have the follwing shape:
* columns:
    * name
    * forward barcode
    * reverse barcorde
```csv
bc1,ATACGATGCTA,GTCGATGTCTGA
b2,GACACACAC,GTCGATTGATG
...
```

## Credit
`nanomux_c` uses `kseq.h` for fastq parsing, and `nob.h`, written by @tsoding (Alexey Kutepov), for overall useful functions!  

## TODO
- [x] trim barcodes
- [x] multi threading
- [] read splitting




