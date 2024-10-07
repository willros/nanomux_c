# nanomux_c

Demultiplex your nanopore (or other) reads! `nanomux` fuzzy matching useful for noisy reads. Written entirely in C, so should compile on most platforms. Depends on `pthreads`. 

In the repo, you can also find `trim` – a small threaded program which you can use to filter out reads with a mean quality and mean length over your decided threshold.

## Quick Start
```bash
$ git clone https://github.com/willros/nanomux_c.git
$ cd nanomux_c/
$ chmod +x build.sh
$ ./build.sh 

# test nanomux
$ ./nanomux -b tests/bc_test.csv -f tests/test.fastq -r 600 -R 2000 -p 200 -k 1 -o new_nanomux -t trim -s split -j 4

# test trim
$ ./trim -i tests/test.fastq -r 2000 -q 20 -t 4
```

## Nanomux
Just run `./nanomux` to get the help message:
```bash
[USAGE]: nanomux 
   -b <barcode>             Path to barcode file.
   -f <fastq>               Path to fastq file.
   -r <read_len_min>        Minimum length of read.
   -R <read_len_max>        Maximum length of read.
   -p <barcode_position>    Position of barcode.
   -k <mismatches>          Number of misatches allowed.
   -o <output>              Name of output folder.
   -t <trim_option>         Trim reads from adapters or not?                      [Options]: trim | notrim.
   -s <split_option>        Split concatenated reads based on adapter occurance?. [Options]: split | nosplit.
   -j <threads>             Number of threads to use.                             Default: 1
```
* `barcode_pos`: where to search for the barcode in the ends. if barcode_pos == 200 and the read length is 1000, the barcodes will be searched for from position 0 -> 200 and 800 -> 1000.
* `k`: allowed number of mismatches
* `trim_option`: To trim the reads to the left and right of the found barcode.
* `split_option`: Try to split reads longer than 2500 if an adapter sequence is found in the middle. One splitted read results in two new reads, with the suffix `_1` and `_2`.

`barcode_file.csv` **MUST** have the follwing shape:
* columns:
    * name
    * forward barcode
    * reverse barcode

It must look like the example below. 
```csv
bc1,ATACGATGCTA,GTCGATGTCTGA
bc2,GACACACAC,GTCGATTGATG
```

## Trim
Run `./trim` to get the help message:
```bash
[USAGE]: trim -i <input> [options]
   -i    <input>             Path of folder or file
   -r    <read_length_min>   Minium length of read.    Default: 1
   -R    <read_length_max>   Minium length of read.    Default: INT_MAX
   -q    <quality>           Minimum quality of read.  Default: 1
   -t    <threads>           Number of threads to use. Default: 1
```

`trim` saves the trimmed reads to a gzipped file with the same name and in the same folder as the original, but with the suffix `.filtered`. 

## Credit
`nanomux_c` uses `kseq.h` for fastq parsing, and `nob.h`, written by @tsoding (Alexey Kutepov), for overall useful functions!  
It also uses `thpool.h` by Johan Hanssen Seferidis.

## TODO
- [x] trim barcodes
- [x] multi threading
- [x] read splitting
- [x] Change to threadpool in nanomux




