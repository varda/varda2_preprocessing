import sys
from os.path import commonprefix
import argparse
import cyvcf2


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


def count_positive_integers(integer_list):
    count = 0
    for integer in integer_list:
        if integer > 0:
            count += 1
    return count


def count_non_ref_alleles(genotypes):
    assert len(genotypes) == 1
    # Assume we have one variant
    first_sample = genotypes[0]
    # The format is [allele_0, allele_1, ... allele_N, <Phasing>]
    without_phasing = first_sample[:-1]

    return count_positive_integers(without_phasing)


def vcf2variants(filename):

    pid_dict = {}
    pid_inc = 1

    # Loop over all variants in a VCF
    for variant in cyvcf2.VCF(filename, gts012=True, lazy=True):

        # shorthands
        ref = variant.REF
        assert len(variant.ALT) == 1
        alt = variant.ALT[0]
        chrom = variant.CHROM

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
        inserted_len = len(inserted)
        if inserted_len == 0:
            inserted = '.'

        # Normalize start and end positions
        norm_start = variant.start + prefix_len
        norm_end = variant.end - suffix_len

        def pgt2phase_set_id(p):
            nonlocal pid_inc
            if p in pid_dict:
                v = pid_dict[p]
            else:
                v = pid_inc
                pid_inc += 1
                pid_dict[p] = v
            return v

        # Determine phase set ID
        # TODO: ideally, if a phase set id only occurs once, because of filtering,
        #       we would like to mark the variant unphased
        pgt_array = variant.format('PGT')
        varda_pid = 0
        if pgt_array:
            if len(pgt_array) != 1:
                print("ERROR: number of phased genotypes is not 1")
                continue

            pgt = pgt_array[0]
            pid = variant.format('PID')[0]

            if pgt == '1|1':
                varda_pid = -1
            elif pgt == '0|1':
                varda_pid = pgt2phase_set_id(f"{chrom}_{pid}_B")
            elif pgt == '1|0':
                varda_pid = pgt2phase_set_id(f"{chrom}_{pid}_A")
            else:
                print('UNEXPECTED PGT: %s' % pgt, file=sys.stderr)

        else:

            # Make sure there isnt more than 1 genotype
            if len(variant.gt_types) != 1:
                print("ERROR: number of gt_types is not 1")
                continue

            # Hetero or Homo?
            if variant.gt_types[0] == 2:
                varda_pid = -1
            else:
                varda_pid = 0

        # Determine the allele_count for the variants
        non_ref_allele_count = count_non_ref_alleles(variant.genotypes)

        print(variant.CHROM, norm_start, norm_end, non_ref_allele_count, varda_pid, inserted_len, inserted, sep='\t')


def main():

    parser = argparse.ArgumentParser(description='Convert VCF to variants.')
    args = parser.parse_args()

    vcf2variants('-')
