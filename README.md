# hicfilters
This README provides documentation for the script `filterHiCReads.py`. The script is intended for flagging HiC reads in the bam file for filtering and QC.


### Author
Chris Nam, at Harvard Medical School.

### Dependencies
Python version >=2.7

```
pip install pysam==0.8.4
pip install biopython
```

### Usage
```
python filterHiCReads.py -i <input.bam> -r <restriction_enzyme_sites.bed.gz> -o <output.filtered.bam> [--verbose]

# If --verbose flag is specified, the script will output intermediate statistics (every 10000 pairs) to stdout.
```

### Output flags
```
YM - whether either read was unmapped or poorly mapped
YL - whether the reads mapped to multiple distal loci. If the reads mapped to > 3 loci, or mapped to 3 loci spaced > 1000 bp apart
YC - whether the reads represent a very close-range contact. If the inter-distance between reads is < 1000 bp
YF - whether the reads mapped to multiple restriction fragments. If the reads mapped to > 3 restriction fragments
YU - whether the reads contain an undigested restriction site. If any of the reads intersected with multiple adjacent restriction fragments
YS - whether the reads mapped to the same restriction fragment
YI - whether the reads come from an insert with anomalous length. If the read pair has an insert size of > 1000 bp
YR - whether either read mapped very far from the nearest restriction site. If either read mapped to > 750 bp of nearest restriction site.
Z1 - 5' position of read 1
Z2 - 5' position of read 2
ZM - position of nearest restriction site downstream of 3' of read 1
ZN - position of nearest restriction site downstream of 3' of read 2
```

