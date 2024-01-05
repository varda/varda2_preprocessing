import pytest
from natsort import natsorted
from vcfphasesets import read_vcf
from vcf2variants import to_varda, trim


@pytest.mark.parametrize("filename, expected", [
    ("tests/data/test_homozygotes.vcf", [("NG_008376.4", 7869, 7870, 3, -1, 1, "T")]),
    ("tests/data/test_unphased.vcf", [("NG_008376.4", 7869, 7870, 2, 0, 1, "T")]),
    ("tests/data/test_triploid.vcf", [
        ("NG_008376.4", 5118, 5119, 1, 2, 1, "T"), ("NG_008376.4", 5156, 5156, 1, 1, 1, "T"),
        ("NG_008376.4", 6754, 6755, 3, -1, 1, "C"), ("NG_008376.4", 7869, 7870, 1, 4, 1, "T"),
        ("NG_008376.4", 7869, 7870, 1, 5, 1, "T"), ("NG_008376.4", 8007, 8008, 1, 3, 1, "A"),
        ("NG_008376.4", 8007, 8008, 1, 4, 1, "A"), ("NG_008376.4", 9199, 9200, 1, 0, 1, "C"),
    ]),
    ("tests/data/test_mixed.vcf", [
        ("NG_008376.4", 7869, 7870, 2, 0, 1, "T"), ("NG_008376.4", 8007, 8008, 3, -1, 1, "A"),
    ]),
    ("tests/data/test_multi_chrom.vcf", [
        ("NG_008376.4", 5118, 5119, 1, 2, 1, "T"), ("NG_008376.4", 5156, 5156, 1, 1, 1, "T"),
        ("NG_008376.X", 5118, 5119, 1, 2, 1, "T"), ("NG_008376.X", 5156, 5156, 1, 1, 1, "T"),
    ]),
])
def test_to_varda(filename, expected):
    _, phase_sets = read_vcf(filename, mod_func=trim)
    assert natsorted(to_varda(phase_sets)) == expected
