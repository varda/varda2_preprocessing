gvcf2coverage
=============

This repository contains two functionally similar implementations of a coverage
extractor from gVCF files.

The Python version is more readable and apt for modification, but the C version
is roughly 16x faster.

They both read the gvcf input (either compressed or uncompressed) from stdin
and output in the following tab seperated format:

`<CHROM> <START> <END> <PLOIDY>`


eg:
```
chrM	0	71	2
chrM	71	72	2
chrM	72	194	2
chrM	194	195	2
chrM	195	301	2
chrM	301	303	2
chrM	303	310	2
chrM	310	311	2
```



Running Python
--------------

python >= 3.6
Virtual environment with cyvcf2.


Building C
----------

Requires htslib.

`make HTSLIB_INCDIR=../../../samtools/include HTSLIB_LIBDIR=../../../samtools/lib`


Running C
---------

`export LD_LIBRARY_PATH=$HTSLIB_LIBDIR`


Recommendation
--------------

Pipe the results of gvcf2coverage(.py) to `bedtools merge` to merge all the
individual adjecent entries.
