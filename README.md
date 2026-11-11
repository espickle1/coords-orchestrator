# coords-orchestrator

Amino acid sequence and coordinate reconciliation tool: one reference frame for protein sequences, per-residue tables and protein structures.

A homologous family has several numberings that look like one integer: the position in a reference
sequence, the row of a per-residue table (one row per residue of *that* sequence), the author residue
number a structure was deposited with, and the column of an alignment. They agree only by coincidence.
Read row 141 of a table for a sequence with a deletion upstream, or ask a structure for "residue 141" when it
was numbered from 501, and you silently read a different residue.

`coords-orchestrator` builds the maps between these numberings **once**, writes them as plain CSV, and gives
you a `Frame` that holds the single reference integer and applies the maps before every read. Where a map is
undefined it returns an explicit *hole* with a reason. It never reads the next row instead.

```
reference position p ──f_i──▶ row of a per-residue table, for any sequence in the family
        │
        └──────g──────▶ (chain, author number, insertion code) of a structure residue
```

Scope: **protein** sequences and **mmCIF** structures. Nucleic-acid chains, ligands and waters in a structure
are read past and dropped.

## Install

```bash
pip install .            # from a checkout; Python >= 3.11
```

Requires numpy, pandas, scipy and biopython. No external aligner is needed.

## Quick start

Describe the project in a `coords.toml` (full reference in [docs/config.md](docs/config.md)):

```toml
[project]
anchor = "A"                    # the group whose references define the frame
out_dir = "maps"

[references.A]                  # one FASTA per group; each record is one protein family
fasta = "refs/A.fasta"
[references.B]
fasta = "refs/B.fasta"

[sequences]                     # optional: sequences you will read per-residue tables for
table = "sequences.csv"         # seq_id, group, family, and sequence | fasta

[structures]                    # optional: mmCIF files
dir = "structures"

[sites]                         # optional: positions given in different numbering schemes
table = "sites.csv"
[sites.schemes]
full_length = 0
mature = 16
```

Build the maps, then read through them:

```bash
python -m coords build coords.toml      # writes maps/
python -m coords info maps              # what was built, and what was flagged
```

```python
from coords import Frame

fr = Frame.load("maps")

fr.rows("kinase", "seq_17", [141, 150])         # 0-based table rows; pd.NA where seq_17 has a hole
fr.struct("kinase", "1abc", 141)                # [ResKey(model, chain, auth_seq_id, ins_code), ...]
fr.ref_pos("1abc", ResKey(1, "A", 49, "A"))     # structure residue -> reference position
fr.neighbors("kinase", "1abc", 141, 8.0,        # residues within 8 Å, each with the row to read
             center="all", seq_ids=["seq_17", "seq_18"])
fr.holes                                        # every hole recorded, with its reason
```

A runnable tour on synthetic data: `python examples/walkthrough.py`.

## The rules

1. **One frame.** Every position that crosses a boundary is a reference position `p` on the anchor
   group's reference of a family. A string position and an author number exist only inside a read. A
   structure residue is a `ResKey`, never a bare integer.
2. **Maps are data.** Each map is built once, written to CSV with the residue letters that justify it, and
   consumed by merge. A reviewer can read a map without reading code.
3. **Holes are explicit and never silently filled.** A lookup returns a value or a status. The default is to
   skip the row and record it in `Frame.holes`; `on_hole="raise"` makes it fatal. There is no pass-through mode.

| Status | Meaning |
|---|---|
| `ok` | A value. |
| `gap_in_sequence` | The sequence has a deletion at that position. |
| `unknown_residue` | The aligned residue is `X` (an unsequenced stretch). A table's row for it describes an unknown token, not residue `p`. |
| `not_observed` | In the structure's SEQRES but has no coordinates (a missing loop). |
| `not_on_structure` | The reference position has no residue in this chain. |
| `other_family` | A neighbour on a chain of another family (a heteromer); its position is in that family's frame. |
| `unregistered_chain` | A neighbour on a chain that matches no reference (an antibody, say). |
| `no_ref_position` | A structure residue the reference has no position for: a tag, a fusion partner, an insertion. |
| `unresolved_scheme` | A site whose numbering scheme is unknown. Excluded, never guessed. |

## What it handles

Indels upstream of the query. Insertion codes (`49` and `49A` are different residues). Offsets that change
inside one chain, so no single `--offset` can fix it. The same region numbered two ways in two structures.
Missing interior loops and chains that do not start at 1. Unobserved residues in both `_pdbx_poly_seq_scheme`
variants. Alternate locations. Homo- and hetero-oligomers. Sequences from a more diverged group, composed
through a reference-to-reference map whose identity is recorded and flagged when it is low. Wildcard `X`
runs. Chains that are not family members at all (an alignment z-score against shuffled copies, not identity
alone, because free end gaps let a short unrelated chain reach a deceptive identity).

## Documentation

| | |
|---|---|
| [docs/alignment_coordinate_systems.md](docs/alignment_coordinate_systems.md) | The problem, stated without the subject matter. |
| [docs/design.md](docs/design.md) | How the maps are made, every output table, the `Frame` API, limits. |
| [docs/config.md](docs/config.md) | The `coords.toml` reference. |
| [docs/existing_tools_survey.md](docs/existing_tools_survey.md) | What already exists, and why this glue layer is separate. |

## Tests

```bash
pip install -e ".[dev]"
pytest
```

The suite runs on synthetic data only (`coords.testing` writes small mmCIF files and a complete toy project
with known ground truth). Each case is one where the right answer differs from the bare index, so code that
passes the integer straight through fails it. `coords.testing.synthetic_frame` builds a `Frame` from sequences
alone, for tests in packages that depend on this one.

## Status

Version 0.1.0. The maps and the `Frame` API are stable enough to use; the CSV `schema_version` is checked on
load and a mismatch asks for a rebuild. The licence has not been chosen yet.
