import argparse
from os.path import commonprefix

from natsort import natsorted

from vcfphasesets import read_vcf


# Strip the prefix from a string
def remove_prefix(text, pref):
    if text.startswith(pref):
        return text[len(pref):]
    return text


# Strip the suffix from a string
def remove_suffix(text, suf):
    if text.endswith(suf):
        return text[:len(text) - len(suf)]
    return text


# Return the common suffix of a list of strings
def common_suffix(entries):
    suf = commonprefix([entry[::-1] for entry in entries])
    return suf[::-1]


def trim(start, end, ref, alt):
    # Find the common prefix of the ref and the alt
    prefix = commonprefix([ref, alt])
    prefix_len = len(prefix)

    # Remove common prefix from both ref and alt
    ref_remain = remove_prefix(ref, prefix)
    alt_remain = remove_prefix(alt, prefix)

    # Now remove common suffix from ref and ald
    suffix = common_suffix([ref_remain, alt_remain])
    suffix_len = len(suffix)

    # Derive inserted string from the remaining alt string
    # If there is no inserted string (len==0), set it to '.'
    inserted = remove_suffix(alt_remain, suffix)

    # Trim start and end positions
    trim_start = start + prefix_len
    trim_end = end - suffix_len

    return trim_start, trim_end, inserted


def to_varda(phase_sets):
    ps_ids = {}
    for (chrom, ps), alleles in phase_sets.items():
        if chrom not in ps_ids:
            ps_ids[chrom] = {(0, 0): 0}

        homozygotes = set.intersection(*(set(allele) for allele in alleles))
        heterozygotes = [set(allele) - homozygotes for allele in alleles]
        unique_variants = set.union(*(set(allele) for allele in heterozygotes))

        for start, end, sequence in homozygotes:
            yield chrom, start, end, len(alleles), -1, len(sequence), sequence or "."

        if len(unique_variants) == 1:
            # We are now a heterozygous unphased variant
            start, end, sequence = unique_variants.pop()
            ploidy = sum(len(allele) for allele in heterozygotes)
            yield chrom, start, end, ploidy, ps_ids[chrom][(0, 0)], len(sequence), sequence or "."
        else:
            for idx, allele in enumerate(heterozygotes):
                for variant in allele:
                    if (ps, idx) not in ps_ids[chrom]:
                        ps_ids[chrom][(ps, idx)] = len(ps_ids[chrom])

                    start, end, sequence = variant
                    yield chrom, start, end, 1, ps_ids[chrom][(ps, idx)], len(sequence), sequence or "."


def main():
    parser = argparse.ArgumentParser(description="Read phase sets from single sample VCF 4.3 file.")
    parser.add_argument("filename", help="VCF file")
    args = parser.parse_args()

    _, phase_sets = read_vcf(args.filename, mod_func=trim)
    for entry in natsorted(to_varda(phase_sets)):
        print(*entry, sep="\t")


if __name__ == "__main__":
    main()
