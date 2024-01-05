import sys
import argparse
from cyvcf2 import VCF


def eprint(*args):
    print(*args, file=sys.stderr)


def gvcf2coverage(threshold, merge, distance):

    vcf = VCF(fname='-', gts012=False, lazy=False, strict_gt=False)

    # eprint(f"samples: {vcf.samples}")
    assert len(vcf.samples) == 1

    # eprint(f"number of seqnames: {len(vcf.seqnames)}")
    assert len(vcf.seqnames) > 0

    first = True

    #
    # Loop over all entries
    #
    for entry in vcf:

        jump = False

        # Depth
        dp = entry.format('MIN_DP')

        if dp is None:
            dp = entry.format('DP')
            if dp is None:
                depth = 0
            else:
                depth = dp[0][0]
        else:
            depth = dp[0][0]

        #
        # If depth is below the threshold, no need to go proceed
        #
        if depth < threshold:
            continue

        #
        # Convenience handles
        #
        chrom = entry.CHROM
        start = entry.start
        end = entry.end
        ploidy = entry.ploidy

        #
        # When we don't merge, just print here and proceed
        #
        if not merge:
            print(chrom, start, end, ploidy, sep="\t")
            continue

        #
        # Open window for first entry
        #
        if first:
            window_start = start
            window_end = end
            window_chrom = chrom
            window_ploidy = ploidy

            first = False
            continue

        #
        # Detect if we need to close and open a new window
        #
        if window_chrom != chrom or window_ploidy != ploidy or window_end + distance < start:
            # Close the window
            print(window_chrom, window_start, window_end, window_ploidy, sep="\t")

            window_start = start
            window_end = end
            window_chrom = chrom
            window_ploidy = ploidy
        else:
            # Extend the window
            window_start = min(window_start, start)
            window_end = max(window_end, end)

    #
    # Always print the last entry when merging, except if there were no entries
    #
    if merge and not first:
        print(window_chrom, window_start, window_end, window_ploidy, sep="\t")


def main():

    parser = argparse.ArgumentParser(description='Convert gVCF to coverage.')
    parser.add_argument('--threshold', '-t', type=int, default=10, help='DP threshold')
    parser.add_argument('--no_merge', '-n', dest='merge', action='store_false', help='Do not merge entries')
    parser.add_argument('--distance', '-d', type=int, default=0, help='Merging distance')
    args = parser.parse_args()

    gvcf2coverage(args.threshold, args.merge, args.distance)
