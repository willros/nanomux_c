# nanomux_c

Demultiplex your nanopore (or other) reads! Implements fuzzy matching useful for noisy reads. Written entirely in C, so should compile on most platforms.


## Quick Start
```bash
$ git clone https://github.com/willros/nanomux_c.git
$ cd nanomux_c/
$ chmod +x build.sh
$ ./build.sh 

# test nanomux
$ ./nanomux tests/bc_test.csv tests/test.fastq 600 2000 200 1 TEST_NANOMUX trim split 4
```

## Information
The fastq parsing still runs in single threaded mode, hence it is quite slow. I will try do implement multithreading for this in the future.

## Help
Just run `./nanomux` to get the help message:
```bash
[USAGE] Add the following arguments in this exact order:
   ./nanomux
   <barcode_file: path>
   <fastq_file: path>
   <read_len_min: int>
   <read_len_max: int>
   <barcode_pos: int>
   <k: int>
   <output_folder: path>
   <trim_option: trim|notrim>
   <split_option: split|nosplit>
   <num_threads: int>
```
* `barcode_pos`: where to search for the barcode in the ends. if barcode_pos == 200 and the read length is 1000, the barcodes will be searched for from position 0 -> 200 and 800 -> 1000.
* `k`: allowed number of mismatches
* `trim_option`: To trim the reads to the left and right of the found barcode.
* `split_option`: Try to split reads longer than 2500 if an adapter sequence is found in the middle. One splitted read results in two new reads, with the suffix `_1` and `_2`.

## Usage
`barcode_file.csv` **MUST** have the follwing shape:
* columns:
    * name
    * forward barcode
    * reverse barcode
```csv
bc1,ATACGATGCTA,GTCGATGTCTGA
bc2,GACACACAC,GTCGATTGATG
...
```

## Credit
`nanomux_c` uses `kseq.h` for fastq parsing, and `nob.h`, written by @tsoding (Alexey Kutepov), for overall useful functions!  
It also uses `thpool.h` by Johan Hanssen Seferidis.

## TODO
- [x] trim barcodes
- [x] multi threading
- [x] read splitting
- [] Change to threadpool in nanomux




