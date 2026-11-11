# Design

This note describes how the maps are made, what a build writes, the `Frame` API, and where the method stops
being trustworthy. [alignment_coordinate_systems.md](alignment_coordinate_systems.md) states the problem
without the tooling; this is the procedure.

## The model

For a protein family there is one **anchor**: the reference sequence of the anchor group. Counting its
residues from 1 gives the **reference position** `p`. Everything else is a partial map into or out of that
frame.

| Map | From | To | Built from |
|---|---|---|---|
| `f_i` | `p` | position `q` in sequence `i` | a pairwise alignment of sequence `i` to the reference of its own group, composed with the group's map to the anchor |
| `g` | `p` | `ResKey(model, chain, auth_seq_id, ins_code)` of a structure residue | the structure's deposited SEQRES (`_pdbx_poly_seq_scheme`) aligned to a reference |
| `h` | `q` | row of a per-residue table | not built: row `q - 1`, valid only for the exact string the table was computed on, which `Frame.check_hash` verifies |

"Partial" is the point. A map is defined where the two sides have a residue and undefined elsewhere; the
undefined places are returned as holes with a reason, never as the neighbouring row.

## How the maps are made

**Alignment.** Biopython `PairwiseAligner`: global, BLOSUM62, gap open -10, extend -0.5, free end gaps (so
fragments align). All of it is configurable and recorded, with the Biopython version, in `manifest.json`,
because a different aligner can place a divergent pair differently. Gaps are read from `Alignment.indices`
(`-1` is a gap), and the aligner is checked not to have altered or dropped a residue.

**`X` is a positional placeholder.** BLOSUM62 scores `X` at -1 against everything, which makes a long run of
`X` (an unsequenced stretch) cheaper to leave unaligned than to place; with free end gaps a terminal run is
dropped entirely and the sequence's flanks are misplaced. `X` is scored `+x_score` (default 1), so the
flanking residues anchor the run and it keeps its length. The row is still a hole (`unknown_residue`): it is the
activation of an unknown token, not of residue `p`. Identity is computed over known residues only.

**Sequences.** Each distinct sequence is aligned once to the reference of its declared `(group, family)`.
Sequences with the same text share one block in `sequence_maps.csv`, keyed by `hash12` (the first 12 hex
characters of the sha256).

**Other groups.** A group other than the anchor gets a reference-to-reference map into the anchor
(`reference_maps.csv`). A sequence or structure of such a group is mapped in two hops: to its own group's
reference, then through that map. The weak alignment is therefore done once, on two full-length references
where it can be inspected, instead of being redone inside every structure. Its identity is recorded; it does
not make a genuinely divergent pair easier to compare.

**Structures.** The SEQRES is aligned, not the observed CA letters, so residues flanking a missing loop are
not misplaced. Both `_pdbx_poly_seq_scheme` variants are read: the 12-column form carries an observed flag
(`pdb_mon_id` is `?` for a residue with no coordinates); the 7-column form does not, so observation comes from
intersecting with `_atom_site`. A chain's family is chosen by alignment, not by the entity description:
the reference with the most identical aligned residues over every (group, family), anchor group first so a tie
goes to the anchor. A chain is *registered* only if it clears three gates: identity, coverage of the chain,
and a z-score of the alignment score against shuffled copies of the chain. The z-score matters: with free end
gaps a short unrelated chain can reach 30 % identity, but its z stays near 0 while a real homologue's is
large, even at 40 % identity. A chain that fails is `unregistered`.

Only the first model is read, one location per atom (blank, else `A`, else the first remaining), and only
polypeptide chains.

## What a build writes

All CSV except the manifest. `ref_pos` is always on the anchor unless a column says otherwise.

| File | One row per | Columns |
|---|---|---|
| `reference/<group>/<family>.fasta` | reference | the designated `S_0` |
| `reference_maps.csv` | non-anchor reference residue | `family, from_group, to_group, from_pos, to_pos, from_aa, to_aa, identical, status` |
| `sequences.csv` | input sequence | `seq_id, group, family, <your extra columns>, fasta, record, seq_hash, length, ref_length, n_aligned, n_unknown, identity, coverage_of_ref, n_inserted, n_gap, quality` (+ `sequence` for inline input) |
| `sequence_maps.csv` | reference position of each distinct sequence | `group, family, hash12, ref_pos, seq_pos, ref_aa, seq_aa, identical, status`; `ref_pos` counts the sequence's own group reference |
| `structures.csv` | structure file | `structure` (file stem), `entry_id`, `path` (relative to the maps), `sha256` |
| `structure_chains.csv` | protein chain | `structure, entry_id, model, chain, asym_id, entity, description, family, group, seqres_len, n_observed, n_mapped, hop_a_identity, anchor_identity, direct_identity, direct_agreement, coverage, score_z, quality, n_ins_codes, scheme_variant` |
| `structure_maps.csv` | SEQRES residue | `structure, model, chain, entity, family, seqres_idx, auth_seq_id, ins_code, ref_pos, struct_aa, ref_aa, identical, observed, status, group_ref_pos, group_ref_aa` |
| `sites_resolved.csv` | site | `family, position, residue, scheme, ref_pos, ref_aa, status, matches_ref` + your extra columns |
| `manifest.json` | | `schema_version`, inputs' sha256, aligner and version, thresholds, reference-to-reference identities, counts, per-chain identities |

**Quality flags.** Sequences: `ok`; `partial` (known aligned residues cover under `min_seq_coverage` of the
reference, so it is trustworthy where defined but mostly holes); `unreliable` (identity under
`min_seq_identity` or fewer than `min_seq_aligned` known residues). Chains: `ok`; `low_identity` (composed
identity to the anchor under `low_chain_identity`); `unregistered` (no reference passed the gates). For a chain of
another group `hop_a_identity` is to its own group's reference, `anchor_identity` is the composed identity
to the anchor, and `direct_identity` / `direct_agreement` compare against aligning the chain straight to the
anchor, as a check on the composition. Flags never delete rows: a flagged chain is still mapped, and it is
yours to decide whether to use it.

## `Frame`

`Frame.load(coords_dir="maps", on_hole="skip")` reads the maps. Positions in and out are always reference
positions; `family` is the protein family; `structure` is a file stem.

| Call | Returns |
|---|---|
| `rows(family, seq_id, ps)` | 0-based row of each position in that sequence's table (`pd.NA` for a hole) |
| `seq_pos(family, seq_id, p)` | the 1-based position, or `None` |
| `row_table(family, seq_id, ps)` | long table: `ref_pos, seq_pos, row, ref_aa, seq_aa, status` |
| `inverse_map(family, seq_id)` | whole-sequence array `string index -> p` (0 = none) |
| `to_anchor(group, family, own_pos)` | a position on another group's reference carried to the anchor |
| `own_to_anchor_array(group, family)` | the whole-reference version of that |
| `struct(family, structure, p, chain=None)` | observed `ResKey`s at `p`, one per copy |
| `ref_pos(structure, reskey)` | reference position of a structure residue |
| `neighbors(family, structure, p, radius, center, copy, seq_ids, cross_family, include_holes)` | residues within `radius` Å, in long format |
| `sites(family=None)` | the resolved site table; unresolved rows are recorded as holes |
| `check_hash(seq_id, sequence= or seq_hash=)` | raises if the string is not the one the map was built from |
| `sequence(seq_id)` | the sequence text, verified against its hash, or `None` |
| `holes`, `clear_holes()`, `record_hole(status, **info)` | the hole log |
| `groups`, `families`, `refs`, `sequences`, `chains`, `smaps` | the loaded tables |

`neighbors` options: `center="CA"` is CA to CA; `"all"` is the smallest heavy-atom distance between two
residues, which finds side-chain contacts a CA-only centre misses. `copy` is `"first"` (default), `"all"` or a
chain id. With `seq_ids`, each neighbour row is crossed with every listed sequence of the neighbour's family
and carries `seq_id, seq_pos, row, seq_aa`. A neighbour on another family's chain (a heteromer) is dropped
as `other_family` by default, because its position is in that family's frame; `cross_family=True` keeps it and
reads it against *its own* family's sequences.

A position outside the reference (`p` not in `1..len`) is a `ValueError`, not a hole: resolve sites through
`Frame.sites()` first. A wrong family for a sequence, or an unknown `seq_id`, is also an error.

## Limits

- **Pairwise, not multiple.** Every sequence is aligned to one reference. A reference that is a poor
  representative of an insertion-rich family will under-map it; choose the anchor accordingly.
- **Divergent groups.** A composed map is only as good as its weakest hop. Check `reference_to_reference` in the
  manifest and `anchor_identity` / `quality` in `structure_chains.csv` before pooling across groups; the map
  says which residue is meant, not that the residues are comparable.
- **Insertion codes** are carried as part of the address everywhere here. Tools downstream that only take an
  integer residue number (some structure-preparation tools mishandle insertion codes) are outside this package.
- **First model only**, one altloc, protein chains only.
- **Aligner.** Biopython's, with the parameters you set. Ties between equally good gap placements are broken
  by the aligner, not by biology; a repeated residue at an indel makes the placement ambiguous.
- **`neighbors`** builds a KD-tree per call and reads the structure lazily once. It is fine for thousands of
  queries, not tuned for millions.
