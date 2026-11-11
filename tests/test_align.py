"""Alignment -> partial map. The maps are the whole content of an alignment; these tests pin that holes
are -1, never a neighbouring index, and that X is a positional placeholder."""
from __future__ import annotations

import numpy as np
import pytest

from coords import align as al
from coords.testing import homolog, random_protein

PRE, POST = random_protein(30, 1), random_protein(30, 2)


def test_identity_alignment_is_the_identity_map():
    s = random_protein(80, 3)
    pa = al.align_pair(s, s)
    np.testing.assert_array_equal(pa.ref_to_q, np.arange(80))
    assert pa.identity == 1.0 and pa.n_aligned == 80


def test_deletion_in_query_is_a_hole_and_shifts_downstream():
    ref = PRE + "AGCDE" + POST
    qry = PRE + "ACDE" + POST                         # G deleted at reference index 31
    pa = al.align_pair(ref, qry)
    assert pa.ref_to_q[31] == -1
    assert pa.ref_to_q[30] == 30 and pa.ref_to_q[32] == 31   # A unmoved, C moved down by one
    assert pa.ref_to_q[-1] == len(qry) - 1


def test_insertion_in_query_has_no_reference_position():
    ref = PRE + "AGCDE" + POST
    qry = PRE + "AGWWCDE" + POST
    pa = al.align_pair(ref, qry)
    assert int((pa.q_to_ref < 0).sum()) == 2
    c = ref.index("C", 30)
    assert pa.ref_to_q[c] == c + 2                                # C moved up by two


def test_alignment_never_alters_residues():
    ref = random_protein(120, 4)
    qry = homolog(ref, 5, 0.3, deletions=[(50, 4)], insertions=[(80, "WWW")])
    pa = al.align_pair(ref, qry)
    for i, j in enumerate(pa.ref_to_q):
        if j >= 0:
            assert pa.q_to_ref[j] == i


def test_unsequenced_run_is_placed_by_flanks_and_keeps_its_length():
    ref = random_protein(100, 6)
    qry = ref[:40] + "X" * 20 + ref[60:]
    pa = al.align_pair(ref, qry)
    np.testing.assert_array_equal(pa.ref_to_q, np.arange(100))     # not left unaligned, not shifted
    assert pa.n_known == 80 and pa.identity == 1.0                 # identity is over known residues only


def test_x_placeholder_is_a_parameter():
    """A terminal X-run is the case free end gaps get wrong: at the stock score the run is dropped."""
    ref = random_protein(100, 6)
    qry = "X" * 30 + ref[40:]                                    # the first 40 residues are unsequenced
    default = al.align_pair(ref, qry)
    stock = al.align_pair(ref, qry, al.AlignerParams(x_score=-1.0))       # the stock BLOSUM62 value
    assert default.ref_to_q[40] == stock.ref_to_q[40] == 30               # the flank anchors it either way
    np.testing.assert_array_equal(default.ref_to_q[10:40], np.arange(30))   # run placed, length kept
    assert (stock.ref_to_q[:40] == -1).all()                              # run left unaligned


def test_params_are_hashable_and_selectable():
    a, b = al.AlignerParams(), al.AlignerParams(matrix="BLOSUM45")
    assert a != b and len({a, b}) == 2
    s = random_protein(60, 7)
    assert al.align_pair(s, s, b).identity == 1.0
    assert "biopython" in al.AlignerParams().as_dict()["backend"]


def test_clean_maps_unknown_characters_to_x():
    assert al.clean("ac-d*e\n") == "ACDE"
    assert al.clean("AC?D") == "ACXD"


def test_empty_inputs_are_all_holes():
    pa = al.align_pair("ACDE", "")
    assert (pa.ref_to_q == -1).all() and pa.n_aligned == 0


def test_compose_and_invert():
    a = np.array([0, 1, -1, 2])          # a -> b
    b = np.array([5, -1, 7])             # b -> c
    np.testing.assert_array_equal(al.compose(a, b), [5, -1, -1, 7])
    np.testing.assert_array_equal(al.invert(a, 4), [0, 1, 3, -1])


def test_composition_equals_the_two_hops():
    ref = random_protein(150, 8)
    mid = homolog(ref, 9, 0.2, deletions=[(30, 2)])
    far = homolog(mid, 10, 0.2, insertions=[(90, "GG")])
    ab = al.align_pair(ref, mid).q_to_ref        # mid -> ref
    bc = al.align_pair(mid, far).q_to_ref        # far -> mid
    far_to_ref = al.compose(bc, ab)
    ok = far_to_ref >= 0
    assert ok.sum() > 100
    for q in np.flatnonzero(ok)[:50]:
        assert far_to_ref[q] == ab[bc[q]]
