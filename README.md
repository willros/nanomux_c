# nanomux_c


## Quick Start
```bash
$ sh build.sh
$ ./nanomux tests/bc_test.csv tests/test.fastq 600 2000 200 1 OUTPUT
```
## Help
Just run `./nanomux` to get the help message:
```bash
usage:
   ./nanomux <barcode_file: path> <fastq_file: path>
   <read_len_min: int> <read_len_max: int>
   <barcode_pos: int> <k: int>
   <output folder: path>
```
* barcode_pos: where to search for the barcode in the ends. 
* k: allowed number of mismatches

## Usage
`barcode_file.csv` **MUST** have the follwing shape:
```csv
#name,forward_bc,reverse_bc
bc1,ATACGATGCTA,GTCGATGTCTGA
b2,GACACACAC,GTCGATTGATG
...
```

## TODO
* trim barcodes
* read splitting




