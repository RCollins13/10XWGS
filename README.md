# 10XWGS
Command-line utilities for processing and analyzing 10X linked-read WGS data

**Contact:** Ryan Collins (rlcollins@g.harvard.edu)

All code copyright (c) 2017 Ryan Collins and is distributed under terms of the MIT license.

## getMolecules.py
Iterates through a 10X bam file and estimates original molecule lengths and coordinates from colinear reads with matching BX tags (i.e. 10X barcodes). 
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
**Usage Notes:**  
1. Input bam file must be coordinate-sorted and indexed.  
2. "Parititioning distance" (option -d / --dist) is the maximum distance permitted between two colinear reads with matching BX tags before considering them to have arisen from independent molecules.
