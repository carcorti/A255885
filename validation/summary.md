# Validation summary

Date: 2026-09-04  
Scope: compact raw GitHub package for OEIS A255885

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
| `paper/A255885_v3.tex` | `748b9d770506e36bff299cb83a3aa702938998c9e3f2e678d70a6d6afc85af48` |

The complete public-file inventory is in `checksums.sha256`.

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

The following checks passed on the compact package assembled on 2026-09-04:

- `sh validation/run.sh quick`: package inventory and checksums, 254-slot data
  reconciliation, 5725 direct witness checks, full bounded production
  regression through `B=2024`, and independent C self-test through `B=2024`.
- `sh validation/run.sh paper`: two settled `pdflatex -draftmode` passes with
  no PDF, error, unresolved reference/citation, or box warning.
- manuscript policy audit: 14/14 gates passed on the exact copied TeX source.
- `cffconvert --validate -i CITATION.cff`: CFF 1.2.0 schema validation passed.
- final-empty-record audit: both public data files end in exactly `LF LF`.
- AddressSanitizer plus UndefinedBehaviorSanitizer: the complete bounded
  production regression and the independent C self-test through `B=2024`
  both passed with leak detection disabled for the OpenMP runtime.
- negative-path checks: removing one terminal LF from a temporary b-file copy
  was rejected; changing a checksummed README byte was rejected; invoking full
  mode without the exact opt-in guard was rejected before creating a run.

```text
QUICK PACKAGE VALIDATION: PASS
PAPER DRAFTMODE CHECK: PASS
FULL FRONTIER REPLAY: NOT RERUN; ARCHIVED PASS RETAINED
```

## Production-code gate at raw-package handoff

```text
whole-trajectory process correctness: PASS
cross-project genericity: PASS
external-v1 private hardening: N/A (metadata-only copyright-name correction)
bounded hardening record and cycle limit: PASS (public wrappers tested locally)
canonical artifact: PASS
explicit user preferences: PASS
publication-readiness: FAIL for final release; PASS for Carlo inspection
manuscript publication-policy audit: PASS
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
review/run freeze readiness: FIX-FIRST before any new official campaign
decision: FIX-FIRST before public release; GO for Carlo inspection
```

Final release remains blocked by the real Zenodo DOI and Carlo's review of this
raw package. The publication source contains one documented metadata-only
correction; no executable statement, result, or paper byte was changed.
