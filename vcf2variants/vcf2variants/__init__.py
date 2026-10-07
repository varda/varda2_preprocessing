import argparse
from os.path import commonprefix

from pysam import VariantFile


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


def read_vcf(filename):
    with VariantFile(filename) as vcf_file:
        # Assume single sample VCF files for now (and forever)
        assert len(vcf_file.header.samples) == 1

        phase_sets = {}
        active_ps_id = 0
        for rec in vcf_file.fetch():
            entry = rec.samples[0]

            ps_id = entry.get("PS")
            if entry.phased:
                ps_id = ps_id or active_ps_id
            if not ps_id:
                ps_id = rec.pos
                active_ps_id = ps_id if entry.phased else 0

            key = rec.rid, ps_id
            if key not in phase_sets:
                phase_sets[key] = [[] for _ in entry["GT"]]

            for idx, gt in enumerate(entry["GT"]):
                if gt:
                    alt = rec.alts[gt - 1]
                    if "<" in alt or ">" in alt:
                        raise ValueError("Cannot deal with symbolic alleles")

                    if alt == "*":
                        continue

                    variant = trim(rec.start, rec.stop, rec.ref, alt)
                    phase_sets[key][idx].append(variant)

        return {key: value for key, value in phase_sets.items() if any(value)}


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
    with VariantFile(args.filename) as vcf:
        id_to_chrom = {value.id: value.name for value in vcf.header.contigs.values()}

    phase_sets = read_vcf(args.filename)

    for entry in sorted(to_varda(phase_sets)):
        ref_id, start, end, ploidy, ps_id, length, sequence = entry
        print(id_to_chrom[ref_id], start, end, ploidy, ps_id, length, sequence, sep="\t")


if __name__ == "__main__":
    main()
