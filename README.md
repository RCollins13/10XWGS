# 10XWGS
Miscellaneous command-line utilities for processing and analyzing 10X linked-read WGS data

## getMolecules.py
Iterates through a 10X bam file and estimates original molecule lengths and coordinates from colinear reads with matching RX SAM tags
```
usage: getMolecules.py [-h] [-d DIST] ibam outfile

Estimate original molecule sizes and coordinates from 10X linked-read WGS
barcodes

positional arguments:
  ibam                  Input bam
  outfile               Output bed file

optional arguments:
  -h, --help            show this help message and exit
  -d DIST, --dist DIST  Molecule partitioning distance in bp (default: 50000)
```
