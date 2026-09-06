# OEIS A255885

This repository accompanies a certified extension of OEIS A255885. For each
integer base `b >= 2`, let `C(b)` be the number of composite integers `c` with
`2 <= c < b` and `b^(c-1) == 1 (mod c^2)`. The sequence term `a(n)` is the
least base for which `C(b) = n`.

This definitive package is release `v1.0.2`. Its public repository is
`https://github.com/carcorti/A255885`, and its persistent archive DOI is
`https://doi.org/10.5281/zenodo.22519077`.

The C17/OpenMP program uses an interval-aware inverse root sieve. The
completed computation covered every base through `B = 2000000000`. A separate
C17 arithmetic engine recomputed the complete domain with different root
construction and counter representations.

## Main result

- The public prefix of 24 terms is extended contiguously through `a(85)`.
- `a(85) = 1244678401`.
- No base through the computed frontier has count 86, so
  `a(86) > 2000000000`.
- Another 23 exact minima beyond the first gap are retained separately.
- Of the 254 inspected indices, 108 have exact minima in the computed domain
  and 146 do not.

## Contents

```text
.gitattributes                   repository text-normalization policy
.gitignore                       generated and local-only exclusions
README.md                        package overview and commands
CITATION.cff                     release and citation metadata
LICENSE                          MIT license
src/a255885.c                    corrected publication source
src/Makefile                     frozen production build and test recipe
data/b255885.txt                 canonical contiguous b-file, n=1..85
data/sparse.txt                  23 exact noncontiguous index/value pairs
results/minima.tsv               normalized public 254-slot minima table
results/results.tsv              byte-exact 254-slot production result
validation/a255885_independent.c frozen independent arithmetic engine
validation/Makefile              independent-engine build recipe
validation/validate_a255885.py   direct bounded Python regression suite
validation/check_data.py         public data/certificate verifier
validation/check_package.py      public inventory and checksum verifier
validation/run.sh                fail-closed public validation runner
validation/config.tsv            byte-exact official run configuration
validation/manifest.tsv          byte-exact official run manifest
validation/certificate.tsv       independent minima and witness certificate
validation/notes.md              evidence scope and size boundary
validation/summary.md            executed checks and frozen identities
validation/checksums.sha256      SHA-256 inventory of every other public file
paper/A255885.tex                monolithic manuscript source
paper/A255885.pdf                compiled manuscript
```

External review files, OEIS editorial material, third-party PDFs, generated
binaries, calibration states, and multi-gigabyte campaign arrays are
intentionally excluded. The large arrays are reproducible outputs, not source
inputs. Their sizes and hashes are preserved in `validation/notes.md` and the
byte-exact configuration/manifest records.

## Requirements

```text
Linux/POSIX environment
GCC-compatible C17 compiler with OpenMP support
GNU make
Python 3.9 or later
standard tools: awk, cmp, df, grep, mktemp, sed, sha256sum
```

TeX Live is needed only for the optional manuscript check.

## Build

From the repository root:

```sh
make -C src release
```

The publication source differs from the reviewed source used for the official
run only in the correction of the copyright-holder name. All executable
statements and the mathematical algorithm are unchanged. The byte-exact run
configuration and manifest retain the historical source hash for provenance.
The `v6` wording retained in the frozen source and bounded validator denotes
their internal review lineage; it is distinct from the repository release
version `v1.0.2` recorded in `CITATION.cff`.

Release `v1.0.2` corrects the package auditor so that the documented public
validation command works both in a normal Git clone and in an extracted
GitHub/Zenodo source archive. The scientific source, build recipe, results,
certificate, data files, manuscript source, and PDF are byte-identical to
release `v1.0.1`.

## Quick public validation

```sh
sh validation/run.sh quick
```

This fail-closed runner verifies the package inventory and hashes; reconciles
the 254 result slots, contiguous and sparse files, run metadata, and
certificate; directly replays all 5725 archived witnesses; runs the bounded
Python definition oracle through `B = 2024`; and runs the independent C
self-test. Generated binaries are removed afterward.

The manuscript can be checked twice without producing a PDF:

```sh
sh validation/run.sh paper
```

## Full frontier reproduction

The full calculation is deliberately not part of the quick check. On the
documented workstation the production invocation took about 916 seconds with
five threads and used about 20.49 GiB peak resident memory; the independent
arithmetic stage took about 403 seconds with 16 threads and used about
16.77 GiB.

The guarded full mode requires an explicit, nonexistent absolute run path and
an opt-in environment variable:

```sh
A255885_ALLOW_FULL=YES sh validation/run.sh full /absolute/path/to/new-run
```

It rebuilds both programs, runs the complete one-segment production sweep,
performs the production verifier, recomputes the domain independently, and
checks the regenerated results and certificate semantics. It retains the
large run directory supplied by the caller; it never overwrites or removes it.

## Production command shape

The official campaign used:

```sh
src/a255885 init RUN_DIR 2000000000 1999999996
OMP_NUM_THREADS=5 OMP_DYNAMIC=FALSE OMP_THREAD_LIMIT=5 src/a255885 next RUN_DIR
src/a255885 status RUN_DIR
src/a255885 verify RUN_DIR
```

`next` performs the complete frontier in this configuration. These commands
document the archived workflow; they are not instructions to extend the
sequence beyond the validated bound.

## Data formats

`data/b255885.txt` contains space-separated `n a(n)` rows with no header,
consecutive indices 1 through 85, and exactly one final empty record.

`data/sparse.txt` uses the same two-column grammar for 23 exact indices beyond
the first gap. It is a support table, not an OEIS b-file, and must not be used
to bridge the unresolved index 86.

`results/results.tsv` preserves the production program's tab-separated
grammar. After its magic line and key/value metadata, the `index base` table
contains all indices 1 through 254; unavailable minima are written as `NA`.
`results/minima.tsv` is the normalized publication view of the same table,
with a comment header beginning with `#` and no production metadata.

`validation/certificate.tsv` contains independent minimum rows (`M`) for all
254 indices and witness rows (`W`) for each of the 84 newly determined
indices. Machine-readable elapsed times retain their original precision.

## Paper and citation

The self-contained manuscript source is `paper/A255885.tex`; its compiled
counterpart is `paper/A255885.pdf`. `CITATION.cff` records release `v1.0.2`,
the public repository, and the persistent Zenodo concept DOI
`10.5281/zenodo.22519077`.

The production source carries the SPDX identifier `MIT`, and the repository
includes the corresponding license file. The source, manuscript, citation
metadata, and license consistently identify Carlo Corti.
