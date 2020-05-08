# Varda2 preprocessing

The Varda2 database stores genomic variants and coverage information. To enable
efficient and meaningful insertion into the database we have defined a set of
preprocessing steps that all participating centers should follow. Note that
these steps are not cast in stone and hopefully will converge as a set of best
practices between the centers.

The information about variants comes from VCF files and the information about
coverage comes from gVCF files. The following sections describe the steps
in more detail.

## Workflow

The following figure depicts the process DAG generated from the Snakemake workflow at https://git.lumc.nl/klinische-genetica/capture-lumc/vcf-to-varda.

![Process Flow](dag.png)

### Variants

To extract variants from the VCF file in a way that Varda can process them, multiple steps are involved:

- trim_alt_and_uncalled:
  - `bcftools view --trim-alt-alleles --exclude-uncalled --output-file {output} {input}`
- split_multi:
  - `bcftools norm --multiallelics - --output {output} {input}`
- exclude_alt_star:
  - `bcftools view --exclude 'ALT==\"*\"' --output-file {output} {input}`

The first step is a pipeline of `bcftools` filtering and normalisation to get
rid of alt-alleles and multi-allelic entries so that we end up with a single
variant per line.


The second step is to take the filtered VCF file and convert it into a Varda
variant file.

- vcf2varda:
  - `vcf2variants < {input} > {output}`

The last step is only required if there are no proper refseq id's used, i.e. only `chr1` or even `1` instead of `NC_000001.10`.

- var_cthreepo:
  - `cthreepo --mapfile h37 --infile {input} --id_from uc --outfile {output} --id_to rs`

This outputs the following tab separated format:
`<CHROM> <START> <END> <PLOIDY> <PHASE SET> <INSERTED LENGTH> <INSERTED SEQUENCE>`

e.g.:
```
NC_000001.10    13656   13658   1       0       0       .
NC_000001.10    13895   13896   1       0       1       A
NC_000001.10    14164   14165   1       0       1       G
NC_000001.10    14672   14673   1       0       1       C
NC_000001.10    14698   14699   1       0       1       G
NC_000001.10    14906   14907   1       0       1       G
```

NB:
- `-1` in `PHASE SET` is homozygous, `0` is unphased
- `.` in `INSERTED SEQUENCE` is no insertion (thus deletion only)


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
