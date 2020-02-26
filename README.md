gvcf2coverage
=============


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
