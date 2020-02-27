gvcf2coverage
=============

This repository contains two functionally similar implementations of a coverage
extractor from gVCF files.

The Python version is more readable and apt for modification, but the C version
is roughly 12x faster.

They both read the gvcf input (either compressed or uncompressed) from stdin
and output on stdout in the following tab separated format:

`<CHROM> <START> <END> <PLOIDY>`

e.g.:
```
chr1    10033   10038   2
chr1    10038   10043   2
chr1    10043   10044   2
chr1    10044   10048   2
chr1    10048   10049   2
chr1    10049   10050   2
chr1    10050   10051   2
chr1    10051   10054   2
```


It is recommended to immediately pipe the results of gvcf2coverage(.py) to
`bedtools merge` to merge all the individual adjecent entries.

Python
------

Requirements

python >= 3.6

Virtual environment with cyvcf2.
e.g.:
```
python3 -m venv venv
source venv/bin/activate
pip install cyvcf2
```



Building C
----------

Requires htslib.

`make HTSLIB_INCDIR=../../../samtools/include HTSLIB_LIBDIR=../../../samtools/lib`


Running C
---------

`export LD_LIBRARY_PATH=$HTSLIB_LIBDIR`


Recommendation
--------------

