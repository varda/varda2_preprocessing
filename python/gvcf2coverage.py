import sys
from cyvcf2 import VCF

threshold = 10

def eprint(*args):
    print(*args, file=sys.stderr)

vcf = VCF(fname='-', gts012=False, lazy=False, strict_gt=False)

eprint(f"samples: {vcf.samples}")
assert len(vcf.samples) == 1

eprint(f"number of seqnames: {len(vcf.seqnames)}")
assert len(vcf.seqnames) > 0

for entry in vcf:
    # Depth
    dp = entry.format('DP')

    if dp is None:
        depth = 0
    else:
        depth = dp[0][0]

    if threshold <= depth:
        print(entry.CHROM, entry.start, entry.end, entry.ploidy, sep="\t")
