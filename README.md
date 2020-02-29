# Varda2 preprocessing

## Coverage

This repository contains two functionally similar implementations of a coverage
extractor from gVCF files.

The Python version is more readable and apt for modification, but the C version
is roughly 12x faster.

They both read the gvcf input (either compressed or uncompressed) from stdin
and output on stdout in the following tab separated format:

`<CHROM> <START> <END> <PLOIDY>`

e.g.:
```
NC_000001.10    10033   10038   2
NC_000001.10    10038   10043   2
NC_000001.10    10043   10044   2
NC_000001.10    10044   10048   2
NC_000001.10    10048   10049   2
NC_000001.10    10049   10050   2
NC_000001.10    10050   10051   2
NC_000001.10    10051   10054   2
```


Both tools by default merge the resulting entries with a default merging
distance of 0.  If merging is disabled, it is recommended to immediately pipe
the results of gvcf2coverage(.py) to `bedtools merge` to merge all the
individual adjecent entries. Note that bedtools will also merge the entries
with a different value in the ploidy column, therefore we opted to do the
merging in the gvcf2coverage tool.

### Python

Requirements

python >= 3.6

Virtual environment with cyvcf2.
e.g.:
```
python3 -m venv venv
source venv/bin/activate
pip install cyvcf2
```

```
usage: gvcf2coverage.py [-h] [--threshold THRESHOLD] [--no_merge] [--distance DISTANCE]
```


### C

Requires HTSLIB library.


On the LUMC Slurm Shark cluster this means:
```
module load library/htslib/1.10.2/gcc-8.3.1
make
```

If your location is non-standard you can pass it like this to the makefile:
```
make HTSLIB_INCDIR=../../../samtools/include HTSLIB_LIBDIR=../../../samtools/lib
```

And for running:
```
export LD_LIBRARY_PATH=$HTSLIB_LIBDIR
```

## Variants

```
bcftools view --trim-alt-alleles --exclude-uncalled <INPUT_VCF> | bcftools norm --multiallelics - | bcftools view --exclude 'ALT=="*"' > <OUTPUT_VCF>
```
