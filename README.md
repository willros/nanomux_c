# nanomux_c


## Quick Start
```bash
$ sh build.sh
$ ./nanomux
```

```bash
usage:
   ./nanomux <barcode_file: path> <fastq_file: path>
   <read_len_min: int> <read_len_max: int>
   <barcode_pos: int> <k: int>
   <output folder: path>
```

### Usage
`barcode_file.csv` **MUST** have the follwing shape:
```csv
#name,forward_bc,reverse_bc
bc1,ATACGATGCTA,GTCGATGTCTGA
b2,GACACACAC,GTCGATTGATG
...
```




