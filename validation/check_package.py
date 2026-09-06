#!/usr/bin/env python3
"""Audit the exact compact GitHub publication tree for A255885."""

from __future__ import annotations

import hashlib
from pathlib import Path


EXPECTED = {
    ".gitattributes",
    ".gitignore",
    "CITATION.cff",
    "LICENSE",
    "README.md",
    "data/b255885.txt",
    "data/sparse.txt",
    "paper/A255885.pdf",
    "paper/A255885.tex",
    "results/minima.tsv",
    "results/results.tsv",
    "src/Makefile",
    "src/a255885.c",
    "validation/Makefile",
    "validation/a255885_independent.c",
    "validation/certificate.tsv",
    "validation/check_data.py",
    "validation/check_package.py",
    "validation/checksums.sha256",
    "validation/config.tsv",
    "validation/manifest.tsv",
    "validation/notes.md",
    "validation/run.sh",
    "validation/summary.md",
    "validation/validate_a255885.py",
}

FROZEN = {
    "src/a255885.c": "7db90f5a737ec5c8deb8c8b03125247e3eb40082b470f4b13d311e980abeaed7",
    "src/Makefile": "5449bb5683720e9992cb89b73da410b2f9c908beea2c0d35d888393ceff0b095",
    "validation/a255885_independent.c": "d7f43d9516179eb1972c598e9b9ff72f8cc00b0f09a887f2284390e7a489f2c4",
    "validation/Makefile": "487545971db58c54b46903bdd697a2a95f2a8663f91b0de1046772d488fc030f",
    "validation/validate_a255885.py": "91d7891754ee7a61b445046970c949707c216c74dcc8d04ca4e7e24054c1ed6d",
    "results/results.tsv": "5553fa5eae302e52a212cf4e8466d7ea4c69dfc975d78537019e50bc59dc18d5",
    "validation/config.tsv": "80b683323184ac042173b5e1779008707dd4e472d3934b62c14d86dcf6e0ab1a",
    "validation/manifest.tsv": "10ed5dbddf5e730c222e221547d7387c35cc7458602130d679a8287a9822823e",
    "data/b255885.txt": "a448a1ca4466c248fc0a201ae35e2f88452f1109ea1eac0ffe007a6448a1b95b",
    "data/sparse.txt": "98e9753b4126b752be8a35d74f48de43313a3163a026129acbb1b261efffb6a4",
    "validation/certificate.tsv": "cb7f7335bc4e88bf523abdb94e370b845ac47735fb9fdca8eb2b4805412bbd33",
    "paper/A255885.pdf": "72afe2a7c666d61b8025a6f49aea847362b0d0f2dc6d070e5dcfc2ac6b912e89",
    "paper/A255885.tex": "80bfd0ffbd8f771c4953bbb8de54769e4128403d483b1b172f2a9212e54bcc1f",
    "results/minima.tsv": "24cb0bc314b666ead3fb7ddaf9ea9daff0afcc521a1939c98a588f658351041a",
}


def fail(message: str) -> None:
    raise AssertionError(message)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_checksums(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text(encoding="ascii").splitlines():
        fields = line.split("  ", 1)
        if len(fields) != 2 or len(fields[0]) != 64 or fields[1] in values:
            fail("malformed or duplicate checksum row")
        int(fields[0], 16)
        values[fields[1]] = fields[0]
    return values


def main() -> int:
    root = Path(__file__).resolve().parent.parent
    actual = {
        path.relative_to(root).as_posix()
        for path in root.rglob("*")
        if path.is_file() and path.relative_to(root).parts[0] != ".git"
    }
    if actual != EXPECTED:
        fail(f"public inventory mismatch: missing={sorted(EXPECTED-actual)} extra={sorted(actual-EXPECTED)}")
    for path in root.rglob("*"):
        if path.relative_to(root).parts[0] == ".git":
            continue
        if path.is_symlink():
            fail(f"symlink not permitted: {path.relative_to(root)}")
        if path.is_file() and len(path.name) > 25:
            fail(f"basename exceeds 25 characters: {path.name}")
        if path.is_file() and path.stat().st_size > 1_000_000:
            fail(f"publication file exceeds 1 MB: {path.relative_to(root)}")
    if sum((root / name).stat().st_size for name in actual) > 1_000_000:
        fail("compact package exceeds 1 MB")
    for name, expected in FROZEN.items():
        if digest(root / name) != expected:
            fail(f"frozen identity mismatch: {name}")

    checksums = read_checksums(root / "validation/checksums.sha256")
    expected_checksum_paths = EXPECTED - {"validation/checksums.sha256"}
    if set(checksums) != expected_checksum_paths:
        fail("checksum inventory coverage mismatch")
    for name, expected in checksums.items():
        if digest(root / name) != expected:
            fail(f"global checksum mismatch: {name}")

    cff = (root / "CITATION.cff").read_text(encoding="utf-8")
    required_cff = (
        "cff-version: 1.2.0",
        "type: software",
        'repository-code: "https://github.com/carcorti/A255885"',
        'doi: "10.5281/zenodo.22519077"',
        'version: "v1.0.2"',
        "license: MIT",
    )
    if any(value not in cff for value in required_cff):
        fail("CITATION.cff required metadata mismatch")
    if cff.count("10.5281/zenodo.22519077") != 4:
        fail("CITATION.cff DOI count mismatch")
    paper = (root / "paper/A255885.tex").read_text(encoding="utf-8")
    if paper.count("10.5281/zenodo.22519077") != 1:
        fail("paper DOI mismatch")
    readme = (root / "README.md").read_text(encoding="utf-8")
    for forbidden in ("\\[", "\\operatorname", "\\varphi", "\\le"):
        if forbidden in readme:
            fail(f"LaTeX-only README token found: {forbidden}")
    forbidden_parts = {"paper_review", "oeis", "__pycache__"}
    for name in actual:
        if forbidden_parts.intersection(Path(name).parts):
            fail(f"private or generated path published: {name}")
    obsolete_doi_sentinel = b"zenodo." + b"x" * 8
    for name in actual:
        if obsolete_doi_sentinel in (root / name).read_bytes():
            fail(f"obsolete Zenodo sentinel remains in final package: {name}")
    print(f"PACKAGE AUDIT PASS: files={len(actual)} bytes={sum((root / name).stat().st_size for name in actual)}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (AssertionError, OSError, ValueError) as error:
        print(f"PACKAGE AUDIT FAIL: {error}")
        raise SystemExit(1)
