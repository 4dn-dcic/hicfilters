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

### Example output
```
SRR1658581.1213 65      chr11   73427756        60      101M    =       122543034       49115279        GATCACTTACACATACCAAACATGACAAAGAAAGCATTTTCCAAAGTATGTCAGCCAGTAAGTCCCAACACATCCAGAAAATCTGTAAAGCCTTTCAATAG   CCCFFFFFHHHGHJJJJIJJIJIIJIHJJJJJJIIIJJJJJJJJJJGIIJJJJIJJIJIIJIGGJJJJJJJJJJHHGHGFFFFFFEEEECEDDDDDDDEE:   NM:i:0  MD:Z:101        AS:i:101
        XS:i:19 YM:i:0  YL:i:0  YC:i:0  YF:i:0  YU:i:0  YS:i:0  YI:i:0  YR:i:0  Z1:i:73427756   Z2:i:122543034  ZM:i:73427895   ZN:i:122543169
SRR1658581.1213 129     chr11   122543034       60      101M    =       73427756        -49115279       ACTTTTTTTTTTTTGAGATGGAGTCTCCCTGTCACCCAGGCTGGAGTCCAGTGGTACAATCTTGGCTCACTGCAACCTTTGCCTCCCGGATTCAAGCGATT   CCCFFFFFHHHHHI<FAABGGCF@=@FHIGGGG@EHEGFFEFCCCCC?@CCC@>@@ACCC>@A@C<@C:AAAC@CCBCAC@ACACCCBB5<>CCDCAB?A8   NM:i:0  MD:Z:101        AS:i:101
        XS:i:48 YM:i:0  YL:i:0  YC:i:0  YF:i:0  YU:i:0  YS:i:0  YI:i:0  YR:i:0  Z1:i:73427756   Z2:i:122543034  ZM:i:73427895   ZN:i:122543169
```

