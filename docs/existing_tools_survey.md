# Has the "alignment coordinate system orchestrator" already been built? Partial pieces exist, but the glue layer does not

**No mature, general-purpose, open-source tool does the whole operation end to end: reference position p → author index a on a structure → radius neighbors (across chains and families) → back to reference positions → rows of every per-sequence table T_i, with explicit holes.** All four layers exist as well-tested building blocks, and a few research tools wire up three of them for their own use cases. The best are MIToS.jl (Julia), EVcouplings' `compare` module (Python) and BioStructMap (Python). None of them gives you a reusable "reference-frame registry" that treats undefined mappings as first-class holes, works on arbitrary predicted models with no SIFTS entry, and handles multi-family complexes. The good news is that the missing glue is small: roughly 150–300 lines of Python on top of Biopython (or Biotite) and gemmi.

## TL;DR
- **Verdict:** People have solved pieces of this many times. Every structural bioinformatician rebuilds the glue privately, which is why "offset/gap/author-number" questions keep appearing on forums. Domain-specific generic numbering schemes (GPCRdb, ANARCI/IMGT, KLIFS) are the closest conceptual match to your orchestrator, but they are hard-wired to one protein family. No general version exists for arbitrary families like phage RBPs and tail fibers.
- **Closest existing tools:**
  - MIToS.jl maps MSA columns to PDB residue numbers via SIFTS, with `""` for missing residues, and computes contacts in MSA-column space.
  - EVcouplings maps alignment indices to PDB indices through SIFTS and computes distance maps.
  - BioStructMap maps reference-aligned data onto a structure and runs radius-based 3D windows.
  - None of them handles your full spec. They rely on SIFTS (useless for AlphaFold/ColabFold models of phage proteins), are tuned for one analysis, or don't expose per-homolog f_i maps with explicit holes.
- **Recommendation:** Build a thin orchestrator on these pieces:
  - Biopython `Bio.Align.Alignment.indices` for f_i (gaps are already encoded as −1).
  - gemmi `align_sequence_to_polymer` (or Biopython `PairwiseAligner`) for the S_0 → (chain, auth_seq_id, ins_code) map.
  - Biopython `NeighborSearch` or Biotite `CellList` for the radius query.
  - pandas nullable dtypes for holes.
  - Use SIFTS/PDBrenum only when the structure is a real PDB entry.

## Key Findings

### Comparison of the best candidates against the four layers

| Tool | Lang / license | Maintenance signal | L1: alignment map f_i | L2: structure numbering map | L3: per-residue table projection | L4: neighbor query + back-translation | Explicit holes? |
|---|---|---|---|---|---|---|---|
| **MIToS.jl** | Julia | v3.11.0 released 7 Oct 2026 (registered in the Julia General registry the same day, PR #170896); ~83 stars | Yes (residue↔column mapping kept even after insert columns are deleted) | Yes, via SIFTS (PDB entries only) | Partial (`msaresidues`) | Partial: `msacontacts` gives contact maps in MSA-column space, not radius queries from a reference position | Yes: `""` for missing PDB residues |
| **EVcouplings (`compare`)** | Python / MIT | v0.2.1, Nov 5 2024; ~249 stars | Yes (`alignment_index_mapping`, `map_indices`) | Yes, via SIFTS + its own remapping | Partial (EC tables) | Partial: distance maps in target numbering, cross-chain via `inter_dists`/`remap_complex_chains` | Partly (NaN distances); tuned to coupling analysis |
| **BioStructMap** | Python / MIT | v0.4.1 on PyPI; 37 stars; low activity | Yes (data must be aligned to a reference) | Yes (BLAST+/pairwise2 alignment of reference to chain) | Yes (data onto structure, B-factor output) | Yes, as 3D sliding-window aggregation (radius) | Not documented; I couldn't verify from source |
| **Biopython** | Python / BSD-style | Very active | Yes: `Alignment.indices` (−1 = gap), `inverse_indices` | Manual (Bio.PDB residue id = (het, resseq, icode)) | No (DIY) | Yes: `NeighborSearch.search(..., level="R")` | Yes for L1 (−1) |
| **Biotite** | Python / BSD | Active (1.6.0 Jan 2026; 1.7.x since) | Yes (alignment trace) | Yes: `use_author_fields` switches auth vs label; `ins_code` annotation | No (DIY) | Yes: `CellList` | DIY |
| **gemmi** | C++/Python / MPL-2.0 | 0.7.5, Mar 2 2026; ~323 stars | Pairwise only | Yes: `align_sequence_to_polymer` (structure-aware gap placement), `label_seq`, `auth_seq_id_to_label` | No | Yes (NeighborSearch) | DIY |
| **SIFTS / PDBe API / PDBrenum** | Data + Python / PDBrenum license unverified | SIFTS updated weekly; PDBrenum DB updated 22 Oct 2026; ~30 stars | No | Yes: residue-level auth ↔ SEQRES ↔ UniProt (PDB entries only) | No | No | PDBrenum pushes unmapped residues to 50000+ (mmCIF) or 5000+ (legacy PDB) and ligands to 60000+; per Faezov & Dunbrack (PLoS One 2021), "Chains with no UniProt in SIFTS are not renumbered" |
| **ChimeraX** | Python-scriptable app | Very active | Yes (Sequence Viewer) | Yes (sequence–structure association) | Yes: header values → `seq_*` residue attributes | Yes (interactive `zone`/selection), not a library API | Yes: unassociated residues get no attribute |
| **ProDy (Evol)** | Python / MIT (per PyPI for 2.6.1) | v2.6.1 (released Aug 19, 2025 on PyPI) | Partial: `alignSequenceToMSA`, `refineMSA` | Partial: `mapOntoChain` | Partial | Contacts possible | DIY |
| **Bio3D** | R / GPL | 2.4-5 (Oct 2024); CRAN reports maintainer email undeliverable | Yes: `pdbaln` gives aligned `resno` matrix (gaps NA) | Yes, author `resno` | Partial (`conserv`) | Yes (distance matrices) | Yes (NA) |
| **dms-viz sitemap** | Python CLI + web / MIT | Small project | No | Yes: `reference_site` → `protein_site` table | Yes (visualization) | No | Sites absent from sitemap aren't shown |
| **GPCRdb / ANARCI / KLIFS** | Web + Python (ANARCI GPLv3) | Mature | Yes, a *generic numbering* registry | Yes (PDB numbering tools) | Yes (residue tables) | Not generic | Yes, by design of the scheme |

### Layer-by-layer assessment

**Layer 1: f_i (reference position ↔ per-sequence position ↔ column).** This is solved in Biopython. `Alignment.indices` "returns a 2D NumPy array with the sequence index of each letter in the alignment, with gaps indicated by -1", and `inverse_indices` maps each letter back to its column.\[1\] Building f_i is a few lines:
1. Take the reference row of `indices`.
2. Keep the columns where it is ≥ 0.
3. Read row i at those columns, treating −1 as a hole.

EVcouplings exposes the same idea as `alignment_index_mapping` / `map_indices`.\[2\] MIToS records each residue's number before deleting insert columns when `generatemapping` is true.\[3\] So you don't need a library for L1; you need to stop reading `T_i[p]` directly.

**Layer 2: reference/sequence position ↔ author index.** This is where people get hurt most. There are two routes, and the right one depends on where your structure came from.

- **Deposited PDB entries.** SIFTS is the authoritative route. It "provides up-to-date residue level mapping between UniProt protein sequences and PDB protein structures", updated weekly, available through the PDBe REST API and embedded in PDBe mmCIF files (`_pdbx_sifts_xref_db`).\[4\]\[5\]
  - PDBrenum rewrites author numbering to UniProt numbering from SIFTS. Its paper (Faezov & Dunbrack, PLoS One 2021) found that 38.6% of the 176,507 PDB entries present on April 26, 2020 contain proteins "not numbered according to a standard for each protein (UniProt)".
  - The MIToS authors justify using SIFTS by noting that aligning PDB sequences with MSA sequences "often yields incorrect alignments in those regions flanking missing residues".\[6\]
  - However, PDBrenum "does NOT process arbitrary models, such as AlphaFold predictions", and SIFTS only covers deposited entries.\[7\]
- **Predicted or custom models.** For AlphaFold/ColabFold models of phage RBPs and tail fibers, and for most phage proteins with no UniProt↔PDB link, you must build the map by alignment.
  - gemmi's `align_sequence_to_polymer` is the best tool for this. It uses inter-atomic distances to put gaps where chain breaks actually are ("Gap positions … are passed to the alignment algorithm as places where the gap can be opened without penalty").\[8\] It also gives you `label_seq` and `auth_seq_id_to_label` conversions.\[9\]
  - Biotite makes the auth/label choice explicit (`use_author_fields`).\[10\] Its maintainers have discussed flipping that default (issue #773), so pin it explicitly.\[11\]

**Layer 3: per-residue data projected onto a shared frame and the structure.** This is done in tool-specific ways:
- **BioStructMap** asks you to supply "sequence-aligned data, a reference sequence and PDB format coordinates". It aligns the reference to each chain (BLAST+ or Bio.pairwise2), maps the data, and writes results to the B-factor column.\[12\]
- **ChimeraX** assigns alignment header values as residue attributes (`seq_conservation`, any `seq_<name>`) "regardless of how they are aligned or numbered (association takes care of the sequence-structure mapping)".\[13\] Residues not associated with the alignment get no value, which is effectively an explicit hole.\[14\]\[15\] This is the most robust *interactive* implementation of your idea, but it is a GUI-centric app rather than a pipeline library.
- **dms-viz** encodes the reference-site → structure-site map as an explicit `sitemap.csv` (`reference_site`, `protein_site`).\[16\] That is a useful *file format* pattern to copy.
- **structuremap** (MannLabs) maps proteomics/PTM sites onto AlphaFold models and runs 3D-proximity analyses.\[17\] It is keyed on UniProt numbering, which works only because AlphaFold DB models are numbered 1..N on the UniProt sequence.

**Layer 4: neighbor query and back-translation.** The geometric part is trivial and mature: Biopython's `NeighborSearch.search(center, radius, level="R")`, Biotite's `CellList`, or gemmi's NeighborSearch.\[18\] The back-translation is not provided by any of them, and that is exactly where your problem lives.
- A Biostars thread ("biopython.PDB residue.get_id() does not match sequence index") shows a user getting neighbor residues back as `resseq` author numbers with "no way to trace" them to sequence indices.\[19\]
- The closest packaged versions of "neighbors in reference space" are BioStructMap's 3D sliding window and MIToS/EVcouplings contact maps indexed by MSA column or target position.\[20\]\[21\] Both are aggregate analyses, not a general query API.

**Domain-specific registries prove the concept.** GPCRdb's structure-based generic numbering places gaps in the alignment so that "residues aligned in sequence are those that align in structure". It also ships lookup tables and a tool that renumbers any user PDB.\[22\]\[23\] ANARCI does the same for antibodies (IMGT/Kabat/Chothia; GPLv3), and KLIFS does it for kinase pockets.\[24\] These are exactly your "single reference-frame coordinate system", but each is hard-coded to one family's HMM or structural template. Nothing equivalent exists for phage receptor-binding proteins.

### Evidence that people keep hitting this exact bug
- On Biostars, users ask how to get "the offset" between PDB numbering (e.g. starting at 669) and data starting at 1. They note gaps and even numbering that "went backward for one part of the chain". The advice given is to align, or keep author numbers and map them yourself.\[25\]
- PDBrenum's paper exists *because* author numbering differs from UniProt in a large share of the PDB. Faezov & Dunbrack report that renumbering changed 49,422 structures (28.0%) "due to differences between UniProt and author residue numbering only", 23,448 (13.3%) because of sequence tags alone, and 18,700 (10.6%) for both reasons.
- Insertion-code bugs recur in published code:
  - In the `prodes` package, the parser ignores PDB column 27, so "CYS A 100 followed by ALA A 100A becomes a single Residue".\[26\]
  - PDBFixer had an issue where insertion codes weren't preserved in `addMissingAtoms`.\[27\]
- MIToS's documentation explicitly warns that naive PDB↔MSA pairwise alignment misplaces residues next to missing segments.\[28\]

Together these show that the shared, robust solution you were hoping for doesn't exist as a general library. Each group re-implements it.

## Recommendations

**1. Don't adopt any one tool wholesale.**
- If your structures were deposited PDB entries of a Pfam family, MIToS.jl (Julia) or EVcouplings' `compare` module (Python) would come closest.
- For ColabFold/AlphaFold models of phage tail fibers, their SIFTS dependence is a dealbreaker.
- BioStructMap is worth reading and borrowing from (radius windows, B-factor export), but it is lightly maintained (v0.4.1, 37 stars).\[29\] I couldn't verify how it handles unaligned residues and insertion codes, so don't trust it for hole semantics without testing.

**2. Build the orchestrator (≈150–300 lines) on these pieces.**
- **Residue key.** A frozen dataclass `ResKey(model, chain, auth_seq_id, ins_code)`. Never use a bare int. This alone eliminates the author-index/reference-position confusion. Read with gemmi or with Biotite `use_author_fields=True` set explicitly.
- **f_i maps (L1).** From Biopython `Alignment.indices`: `ref_cols = indices[0] >= 0`, then `f_i = indices[i][ref_cols]` with −1 → `pd.NA`. Store these as a `pandas` DataFrame indexed by reference position p, one nullable `Int64` column per sequence.
- **Structure map (L2), one per (structure, chain) → (family, sequence id).**
  - Align S_0 (or, better, the family member the structure was built from, then compose through that member's f) to the chain's observed residues with `gemmi.align_sequence_to_polymer`.
  - Emit `p → ResKey` with holes where the chain lacks the residue, and store the inverse `ResKey → p`.
  - For deposited entries, cross-check against SIFTS (`_pdbx_sifts_xref_db` in PDBe mmCIF, or the PDBe API).
- **Multi-chain, multi-family complexes.** Keep a registry `chain → (family alignment, row id)`. Homotrimeric tail fibers map three chains to the same row. A partner protein (e.g. a chaperone or baseplate component) maps to its own alignment. Neighbors on unregistered chains come back as `ResKey` with `p = NA`, never silently dropped.
- **Query (L4).** `query(p, radius)`:
  1. `ResKey = struct_map[p]` (return an empty result with an explicit reason if undefined).
  2. Run `NeighborSearch` (Biopython) or `CellList` (Biotite) at residue level. Decide explicitly whether the center is the CA or all atoms of the residue, since a CA-only center misses some contacts.\[30\]
  3. Map each neighbor back to `(family, p_j)` via the inverse map.
  4. Return a long-format table with one row per neighbor per sequence i: `ResKey`, `p_j`, `f_i(p_j)` (nullable), `T_i` value (NA if unmapped), and `E` value looked up by `ResKey`. Include a `status` column (`ok` / `not_on_structure` / `gap_in_S_i` / `other_family`), following the explicit `mapping_status` pattern used by proptm3d.\[31\]
- **Electrostatics/surface table E.** Key APBS/PDB2PQR or SASA outputs on `ResKey`, not on row number. PDB2PQR/PQR round-trips can renumber or drop insertion codes, so validate that the keys round-trip.
- **Tests.** Write property tests on synthetic alignments with indels before and after p. Add a structure with an insertion code, a missing loop, and a non-1 start. Each of these is a documented failure mode above.

**3. Optional conveniences.** Export per-reference-position values as a ChimeraX `defattr` file or as Jalview/ChimeraX alignment headers, so ChimeraX's association logic independently cross-checks your mapping visually.\[13\] Use a dms-viz-style `sitemap.csv` as the human-readable, versioned record of the p ↔ ResKey map.\[32\]

## Caveats
- I searched broadly but not exhaustively; GitLab/Zenodo/Bioconductor coverage was thin. A small, unpublicized repository implementing a general "residue registry" could exist. I found none with meaningful adoption.
- Version, star and date figures come from PyPI/GitHub/CRAN snapshots and some secondary aggregators (libraries.io for EVcouplings and gemmi stars). Treat them as approximate.
- Not verified:
  - BioStructMap's handling of gaps and insertion codes. Bio.PDB residue IDs include the insertion code, so it *may* preserve them, but I didn't read its source.
  - PDBrenum's licence. ProDy 2.6.1 is listed as MIT on PyPI, but its 1.x releases were GPLv3, so check the licence of whichever version you pin.
  - Release dates are now confirmed: ProDy 2.6.1 on Aug 19, 2025 (PyPI), and MIToS v3.11.0 on 7 Oct 2026, registered in the Julia General registry the same day (PR #170896).
- Biotite's auth-vs-label default has been under discussion (issue #773), so always pass `use_author_fields` explicitly.\[11\]
- Bio3D is still on CRAN, but CRAN flags its maintainer address as undeliverable, which is a maintenance risk for an R route.\[33\]
- The alignment-based structure map (L2 for predicted models) is only as good as the S_0-to-chain alignment. For divergent tail-fiber homologs, align the structure's *own* source sequence into the MSA (it should already be an MSA row) and compose maps, rather than aligning S_0 directly to the structure.

## Sources

1. [Sequence alignments — Biopython 1.85 documentation](https://biopython.org/docs/1.85/Tutorial/chapter_align.html)
2. [EVcouplings Documentation Release 0.0.1 Thomas Hopf Nov 06, 2020](https://app.readthedocs.org/projects/evcouplings/downloads/pdf/latest/)
3. [MSA · MIToS](https://docs.juliahub.com/MIToS/PuGh5/2.10.0/MSA/)
4. [Unified access to up-to-date residue-level annotations from UniProtKB and other biological databases for PDB data](https://www.nature.com/articles/s41597-023-02101-6)
5. [Unified access to up-to-date residue-level annotations from UniProt and other biological databases for PDB data via PDBx/mmCIF files](https://www.biorxiv.org/content/10.1101/2022.08.10.503473.full.pdf)
6. [MIToS.jl: mutual information tools for protein sequence analysis in the Julia language](https://academic.oup.com/bioinformatics/article/33/4/564/2608634)
7. [Renumbering PDB files - Proteopedia, life in 3D](https://proteopedia.org/wiki/index.php/Renumbering_PDB_files)
8. [Structure analysis - Gemmi 0.7.5 documentation](https://gemmi.readthedocs.io/en/stable/analysis.html)
9. [Gemmi 0.7.6-dev documentation](https://gemmi.readthedocs.io/)
10. [get\_structure — Biotite](https://www.biotite-python.org/1.5.0/apidoc/biotite.structure.io.pdbx.get_structure.html)
11. [Make \`use\_author\_fields\` default to \`True\` · Issue #773 · biotite-dev/biotite · GitHub](https://github.com/biotite-dev/biotite/issues/773)
12. [BioStructMap: a Python tool for integration of protein structure and sequence-based features](https://academic.oup.com/bioinformatics/article/34/22/3942/5042171)
13. [\[chimerax-users\] User defined color values for each residue/position](https://www.rbvi.ucsf.edu/pipermail/chimerax-users/2020-November/001747.html)
14. [Showing Sequence Conservation in Chimera](https://www.cgl.ucsf.edu/chimera/data/tutorials/systems/outline.html)
15. [Multalign Viewer (Sequence)](https://www.cgl.ucsf.edu/chimera/docs/ContributedSoftware/multalignviewer/multalignviewer.html)
16. [Data Requirements](https://dms-viz.github.io/dms-viz-docs/preparing-data/data-requirements/)
17. [GitHub - MannLabs/structuremap: Python package for investigating the structural context of PTMs · GitHub](https://github.com/MannLabs/structuremap)
18. [Bio.PDB.NeighborSearch module — Biopython 1.88.dev0 documentation](https://biopython.org/docs/dev/api/Bio.PDB.NeighborSearch.html)
19. [biopython.PDB residue.get\_id() does not match sequence index](https://www.biostars.org/p/401105/)
20. [Proteome-wide mapping of immune features onto Plasmodium protein three-dimensional structures - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC5847524/)
21. [Proteome-wide mapping of immune features onto Plasmodium protein three-dimensional structures](https://www.nature.com/articles/s41598-018-22592-3)
22. [GPCRdb: an information system for G protein-coupled receptors](https://academic.oup.com/nar/article/44/D1/D356/2502664)
23. [GPCRdb: the G protein‐coupled receptor database](https://bpspubs.onlinelibrary.wiley.com/doi/10.1111/bph.13509)
24. [ANARCI: antigen receptor numbering and receptor classification](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4708101/)
25. [Pdb Coordinate Question](https://www.biostars.org/p/5241/)
26. [Insertion Codes Dropped · Issue #14 · datacatalysis/prodes · GitHub](https://github.com/datacatalysis/prodes/issues/14)
27. [Insertion codes not preserved in addMissingAtoms · Issue #187 · openmm/pdbfixer · GitHub](https://github.com/openmm/pdbfixer/issues/187)
28. [MIToS.jl/docs/src/SIFTS.md at master · diegozea/MIToS.jl · GitHub](https://github.com/diegozea/MIToS.jl/blob/master/docs/src/SIFTS.md)
29. [GitHub - andrewguy/biostructmap: A Python package for mapping sequence aligned data onto protein structures](https://github.com/andrewguy/biostructmap)
30. [Adopt upstream #5380: NeighborSearch residue examples in the PDB tutorial by dbolser · Pull Request #321 · BioAGI-org/biop-ai · GitHub](https://github.com/BioAGI-org/biop-ai/pull/321)
31. [GitHub - prolfqua/proptm3d: visualize ptm's directly on the 3d structure along with log2FC · GitHub](https://github.com/prolfqua/proptm3d)
32. [pypi.org](https://pypi.org/project/configure-dms-viz/0.1.3/)
33. [CRAN: Package bio3d - R Project](https://cran.r-project.org/package=bio3d)
