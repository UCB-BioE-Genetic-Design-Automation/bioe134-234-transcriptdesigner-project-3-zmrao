import pytest
from genedesign.checkers.internal_rbs_checker import InternalRBSChecker


@pytest.fixture
def checker():
    c = InternalRBSChecker()
    c.initiate()
    return c


def test_internal_rbs_detected(checker):
    # SD motif + spacer + start codon should be flagged.
    seq = "AAAGGAGGAAAAATGAAACCC"
    passed, hit = checker.run(seq)
    assert passed is False
    assert hit is not None
    assert "ATG" in hit


def test_internal_rbs_not_detected_without_start(checker):
    # SD motif present but no start codon in spacer window.
    seq = "AAAGGAGGAAAAAACCCGGG"
    passed, hit = checker.run(seq)
    assert passed is True
    assert hit is None


def test_internal_rbs_not_detected_random(checker):
    seq = "AACCTTGGCCAATTCGATCGG"
    passed, hit = checker.run(seq)
    assert passed is True
    assert hit is None
