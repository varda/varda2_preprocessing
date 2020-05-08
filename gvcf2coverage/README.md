# gvcf2coverage 

## Requirements
Requires HTSLIB library.

## Installation
On the LUMC Slurm Shark cluster this means:
```
module load library/htslib/1.10.2/gcc-8.3.1
make
```

If your location is non-standard you can pass it like this to the makefile:
```
make HTSLIB_INCDIR=../../../samtools/include HTSLIB_LIBDIR=../../../samtools/lib
```

## Usage
And for running:
```
export LD_LIBRARY_PATH=$HTSLIB_LIBDIR
```
