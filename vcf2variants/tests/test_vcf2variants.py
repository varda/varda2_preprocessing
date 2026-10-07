import pytest
from vcf2variants import read_vcf,to_varda


@pytest.mark.parametrize("filename, expected", [
    ("tests/data/test_homozygotes.vcf", [
        (0, 7869, 7870, 3, -1, 1, "T")
    ]),
    ("tests/data/test_unphased.vcf", [
        (0, 7869, 7870, 2, 0, 1, "T")
    ]),
    ("tests/data/test_triploid.vcf", [
        (0, 5118, 5119, 1, 2, 1, "T"),
        (0, 5156, 5156, 1, 1, 1, "T"),
        (0, 6754, 6755, 3, -1, 1, "C"),
        (0, 7869, 7870, 1, 4, 1, "T"),
        (0, 7869, 7870, 1, 5, 1, "T"),
        (0, 8007, 8008, 1, 3, 1, "A"),
        (0, 8007, 8008, 1, 4, 1, "A"),
        (0, 9199, 9200, 1, 0, 1, "C"),
    ]),
    ("tests/data/test_mixed.vcf", [
        (0, 7869, 7870, 2, 0, 1, "T"),
        (0, 8007, 8008, 3, -1, 1, "A"),
    ]),
    ("tests/data/test_multi_chrom.vcf", [
        (0, 5118, 5119, 1, 2, 1, "T"),
        (0, 5156, 5156, 1, 1, 1, "T"),
        (1, 5118, 5119, 1, 2, 1, "T"),
        (1, 5156, 5156, 1, 1, 1, "T"),
    ]),
])
def test_to_varda(filename, expected):
    phase_sets = read_vcf(filename)
    assert sorted(to_varda(phase_sets)) == expected
