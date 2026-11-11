"""Pairwise alignment -> partial map ``f : reference position -> query position``.

Backend: Biopython ``PairwiseAligner`` (global, a BLOSUM/PAM matrix, affine gaps, free end gaps by
default so partial sequences are tolerated). ``Alignment.indices`` gives the gap encoding directly
(``-1`` = gap), so a hole is a ``-1`` in the output, never a fabricated row.

``X`` is a positional placeholder. BLOSUM matrices score X at about -1 against everything, which
makes a long run of X (an unsequenced stretch) cheaper to leave unaligned than to place. Scoring X at
``+x_score`` lets the flanking residues anchor the run so it keeps its length.

Every map built from these alignments records the parameters and the Biopython version in
``manifest.json``, because two aligners can disagree on divergent pairs.
"""
from __future__ import annotations

import re
from dataclasses import asdict, dataclass
from typing import Dict

import numpy as np


@dataclass(frozen=True)
class AlignerParams:
    """Settings of the global aligner. Frozen so it can key a cache."""

    matrix: str = "BLOSUM62"          # any name Bio.Align.substitution_matrices.load accepts
    open_gap: float = -10.0
    extend_gap: float = -0.5
    end_gap: float = 0.0              # 0 = free end gaps (partial sequences are tolerated)
    x_score: float = 1.0              # score of X against anything (positional placeholder)

    def as_dict(self) -> Dict[str, object]:
        d = asdict(self)
        d["backend"] = "biopython.PairwiseAligner"
        d["mode"] = "global"
        return d


DEFAULT_PARAMS = AlignerParams()
_CACHE: Dict[AlignerParams, object] = {}


def _aligner(params: AlignerParams):
    if params not in _CACHE:
        from Bio import Align
        from Bio.Align import substitution_matrices

        a = Align.PairwiseAligner()
        a.mode = "global"
        m = substitution_matrices.load(params.matrix)
        alpha = m.alphabet
        arr = np.array(m, dtype=float)
        if "X" in alpha:
            xi = alpha.index("X")
            arr[xi, :] = params.x_score
            arr[:, xi] = params.x_score
        a.substitution_matrix = substitution_matrices.Array(alphabet=alpha, dims=2, data=arr)
        a.open_gap_score = params.open_gap
        a.extend_gap_score = params.extend_gap
        a.end_gap_score = params.end_gap
        _CACHE[params] = a
    return _CACHE[params]


def biopython_version() -> str:
    import Bio

    return Bio.__version__


def clean(seq: str) -> str:
    """Upper-case, strip whitespace/gaps/stop. Anything outside the amino-acid alphabet becomes X."""
    s = re.sub(r"[\s\-*]", "", str(seq).upper())
    return re.sub(r"[^ACDEFGHIKLMNPQRSTVWYBZXJUO]", "X", s)


@dataclass
class PairAlignment:
    """Alignment of ``query`` against ``ref``.

    ``ref_to_q`` has one entry per reference residue: the 0-based query index, or ``-1`` when the
    query has a gap there. ``q_to_ref`` is the inverse (``-1`` = insertion relative to the ref).
    """

    ref_to_q: np.ndarray
    q_to_ref: np.ndarray
    score: float
    n_aligned: int
    n_identical: int
    n_known: int = -1         # aligned columns where neither side is X; identity is taken over these

    @property
    def identity(self) -> float:
        n = self.n_aligned if self.n_known < 0 else self.n_known
        return self.n_identical / n if n else float("nan")

    def coverage_of_ref(self, n_ref: int) -> float:
        return self.n_aligned / n_ref if n_ref else float("nan")

    def coverage_of_query(self, n_q: int) -> float:
        return self.n_aligned / n_q if n_q else float("nan")


def align_pair(ref: str, query: str, params: AlignerParams = DEFAULT_PARAMS) -> PairAlignment:
    """Global alignment of ``query`` to ``ref``. Raises if the aligner altered or dropped a residue."""
    ref, query = clean(ref), clean(query)
    if not ref or not query:
        return PairAlignment(np.full(len(ref), -1), np.full(len(query), -1), float("-inf"), 0, 0, 0)
    aln = _aligner(params).align(ref, query)[0]
    idx = np.asarray(aln.indices)
    r, q = idx[0], idx[1]
    both = (r >= 0) & (q >= 0)
    ref_to_q = np.full(len(ref), -1, dtype=np.int64)
    q_to_ref = np.full(len(query), -1, dtype=np.int64)
    ref_to_q[r[both]] = q[both]
    q_to_ref[q[both]] = r[both]
    if int((r >= 0).sum()) != len(ref) or int((q >= 0).sum()) != len(query):
        raise RuntimeError("alignment dropped residues")
    ident = int(sum(ref[i] == query[j] for i, j in zip(r[both], q[both]) if ref[i] != "X"))
    known = int(sum(ref[i] != "X" and query[j] != "X" for i, j in zip(r[both], q[both])))
    return PairAlignment(ref_to_q, q_to_ref, float(aln.score), int(both.sum()), ident, known)


def compose(first: np.ndarray, second: np.ndarray) -> np.ndarray:
    """Compose ``first: a -> b`` with ``second: b -> c`` (arrays of 0-based indices, -1 = hole)."""
    out = np.full(len(first), -1, dtype=np.int64)
    ok = first >= 0
    out[ok] = np.where(first[ok] < len(second), second[np.clip(first[ok], 0, max(len(second) - 1, 0))], -1)
    return out


def invert(m: np.ndarray, n_target: int) -> np.ndarray:
    """Invert a partial injective map ``a -> b`` into ``b -> a`` (-1 = no preimage)."""
    out = np.full(n_target, -1, dtype=np.int64)
    ok = m >= 0
    out[m[ok]] = np.flatnonzero(ok)
    return out
