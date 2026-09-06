# Validation summary

Date: 2026-09-06  
Scope: compact corrective GitHub/Zenodo release package `v1.0.2` for OEIS A255885

Release `v1.0.2` changes only packaging validation and release metadata. It
removes an erroneous dependency on the checkout-directory basename, excludes
the root `.git/` metadata directory from the public-file inventory, and adds
the tracked `.gitattributes` file to that inventory. The scientific source,
build recipe, results, certificate, data files, manuscript source, and PDF are
byte-identical to release `v1.0.1`.

## Frozen identities

| Artifact | SHA-256 |
|---|---|
| `src/a255885.c` | `7db90f5a737ec5c8deb8c8b03125247e3eb40082b470f4b13d311e980abeaed7` |
| `src/Makefile` | `5449bb5683720e9992cb89b73da410b2f9c908beea2c0d35d888393ceff0b095` |
| `validation/a255885_independent.c` | `d7f43d9516179eb1972c598e9b9ff72f8cc00b0f09a887f2284390e7a489f2c4` |
| `validation/Makefile` | `487545971db58c54b46903bdd697a2a95f2a8663f91b0de1046772d488fc030f` |
| `validation/validate_a255885.py` | `91d7891754ee7a61b445046970c949707c216c74dcc8d04ca4e7e24054c1ed6d` |
| `results/results.tsv` | `5553fa5eae302e52a212cf4e8466d7ea4c69dfc975d78537019e50bc59dc18d5` |
| `results/minima.tsv` | `24cb0bc314b666ead3fb7ddaf9ea9daff0afcc521a1939c98a588f658351041a` |
| `data/b255885.txt` | `a448a1ca4466c248fc0a201ae35e2f88452f1109ea1eac0ffe007a6448a1b95b` |
| `data/sparse.txt` | `98e9753b4126b752be8a35d74f48de43313a3163a026129acbb1b261efffb6a4` |
| `validation/certificate.tsv` | `cb7f7335bc4e88bf523abdb94e370b845ac47735fb9fdca8eb2b4805412bbd33` |
| `paper/A255885.tex` | `80bfd0ffbd8f771c4953bbb8de54769e4128403d483b1b172f2a9212e54bcc1f` |
| `paper/A255885.pdf` | `72afe2a7c666d61b8025a6f49aea847362b0d0f2dc6d070e5dcfc2ac6b912e89` |

The checksum manifest covers every other public file; it excludes itself.

The official run manifest retains the pre-correction source hash
`eab592a8c6a3caf75fa868724765800347516b23b65cc843c6d0a82c205aa1a5`.
The publication source differs only in the corrected copyright-holder name;
all executable statements are unchanged.

## Expected numerical summary

- Bound: `B=2000000000`.
- Complete modulus frontier: `1999999999`.
- Result slots: 254.
- Populated slots: 108; empty slots: 146.
- New exact indices: 84.
- Contiguous endpoint: 85.
- First unresolved index: 86, with `a(86)>2000000000`.
- Additional exact sparse indices: 23.
- Direct archived witness checks: 5725.
- Last contiguous row: `85 1244678401`.

## Public commands

```sh
sh validation/run.sh quick
sh validation/run.sh paper
```

The guarded full-frontier command is documented in the root README. It was
not rerun during compact packaging because that would be a new full campaign;
the archived full validation status is `FINAL ARITHMETIC VALIDATION PASS`.

## Packaging verification state

The following checks passed on the final package assembled on 2026-09-06:

- `sh validation/run.sh quick`: package inventory and checksums, 254-slot data
  reconciliation, 5725 direct witness checks, full bounded production
  regression through `B=2024`, and independent C self-test through `B=2024`.
- `sh validation/run.sh paper`: two settled `pdflatex -draftmode` passes over
  the final unversioned source with no error, unresolved reference/citation,
  or box warning.
- release PDF: two settled normal `pdflatex` passes produced the 10-page
  `paper/A255885.pdf`; the final log and rendered-page inspection found no
  error, warning, unresolved reference/citation, clipped content, overlapping
  content, or split table.
- manuscript policy audit: 14/14 gates passed on the exact copied TeX source.
- `cffconvert --validate -i CITATION.cff`: CFF 1.2.0 schema validation passed.
- distribution-context matrix: the package audit passed from an arbitrarily
  named directory and from a clone-like tree containing root `.git/` metadata;
  the complete quick-validation command passed after ZIP creation and
  extraction under a GitHub-style source-archive directory name.
- final-empty-record audit: both public data files end in exactly `LF LF`.
- AddressSanitizer plus UndefinedBehaviorSanitizer: the complete bounded
  production regression and the independent C self-test through `B=2024`
  both passed with leak detection disabled for the OpenMP runtime.
- negative-path checks: removing one terminal LF from a temporary b-file copy
  was rejected; changing a checksummed README byte was rejected; invoking full
  mode without the exact opt-in guard was rejected before creating a run.

```text
QUICK PACKAGE VALIDATION: PASS
PAPER SOURCE/PDF CHECK: PASS
FULL FRONTIER REPLAY: NOT RERUN; ARCHIVED PASS RETAINED
```

## Production-code gate at final-package handoff

```text
whole-trajectory process correctness: PASS
cross-project genericity: PASS
external-v1 private hardening: N/A (no executable scientific-code change)
bounded hardening record and cycle limit: PASS (public wrappers tested locally)
canonical artifact: PASS
explicit user preferences: PASS
publication-readiness: PASS
manuscript publication-policy audit: PASS
personal identity integrity: PASS (Carlo Corti)
validation identity stability: PASS (run and publication hashes separated)
five-boundary separation: PASS
sequence-agnostic endpoint discovery: PASS
structural/numeric capacity: PASS
measured hardware-envelope projection: PASS
operational segments <=90 minutes: PASS
campaign artifact/filename budget: PASS
external-review attachment budget: N/A
external-review version lineage: PASS with documented metadata-only correction
independent-validator reproducibility: PASS, full replay available but not rerun
qualified OpenMP TSan controls: inherited; not rerun during packaging
review/run freeze readiness: PASS for publication; any new campaign requires a
separate authorization and gate
decision: GO for publication
```

The final metadata consistently records the public repository, the persistent
Zenodo concept DOI `10.5281/zenodo.22519077`, and release `v1.0.2`. The
scientific code, results, b-file, sparse data, run metadata, and certificate
remain unchanged from the validated package.
