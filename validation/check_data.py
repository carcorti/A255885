#!/usr/bin/env python3
"""Fail-closed checks for the compact public A255885 data package."""

from __future__ import annotations

import argparse
import hashlib
import math
from pathlib import Path


BOUND = 2_000_000_000
INDEX_CAP = 254
CONTIGUOUS_END = 85
EXPECTED_NON_NA = 108
EXPECTED_SPARSE = 23
EXPECTED_NEW = 84
EXPECTED_WITNESSES = 5_725


def fail(message: str) -> None:
    raise AssertionError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require_final_blank(path: Path) -> None:
    data = path.read_bytes()
    if not data.endswith(b"\n\n") or data.endswith(b"\n\n\n"):
        fail(f"{path.name} must end in exactly two LF bytes")


def parse_pairs(path: Path) -> list[tuple[int, int]]:
    require_final_blank(path)
    rows: list[tuple[int, int]] = []
    for lineno, line in enumerate(path.read_text(encoding="ascii").splitlines(), 1):
        if not line:
            continue
        fields = line.split()
        if len(fields) != 2 or not all(field.isdecimal() for field in fields):
            fail(f"malformed pair at {path.name}:{lineno}")
        rows.append((int(fields[0]), int(fields[1])))
    return rows


def parse_results(path: Path) -> tuple[dict[str, str], list[int | None]]:
    lines = path.read_text(encoding="ascii").splitlines()
    if not lines or lines[0] != "A255885-RESULTS\t1":
        fail("results magic/version mismatch")
    try:
        marker = lines.index("index\tbase")
    except ValueError:
        fail("results table marker missing")
    header: dict[str, str] = {}
    for line in lines[1:marker]:
        fields = line.split("\t")
        if len(fields) != 2 or fields[0] in header:
            fail("malformed or duplicate results metadata")
        header[fields[0]] = fields[1]
    rows: list[int | None] = []
    for expected, line in enumerate(lines[marker + 1 :], 1):
        fields = line.split("\t")
        if len(fields) != 2 or fields[0] != str(expected):
            fail("results index order mismatch")
        rows.append(None if fields[1] == "NA" else int(fields[1]))
    if len(rows) != INDEX_CAP:
        fail(f"expected {INDEX_CAP} result rows, got {len(rows)}")
    return header, rows


def parse_public_minima(path: Path) -> list[int | None]:
    require_final_blank(path)
    lines = path.read_text(encoding="ascii").splitlines()
    if not lines or lines[0] != "# index\tbase":
        fail("normalized minima header mismatch")
    rows: list[int | None] = []
    for expected, line in enumerate((line for line in lines[1:] if line), 1):
        fields = line.split("\t")
        if len(fields) != 2 or fields[0] != str(expected):
            fail("normalized minima index order mismatch")
        rows.append(None if fields[1] == "NA" else int(fields[1]))
    if len(rows) != INDEX_CAP:
        fail("normalized minima row-count mismatch")
    return rows


def parse_key_values(path: Path, magic: str) -> dict[str, str]:
    lines = path.read_text(encoding="ascii").splitlines()
    if not lines or lines[0] != magic:
        fail(f"{path.name} magic/version mismatch")
    values: dict[str, str] = {}
    for line in lines[1:]:
        fields = line.split("\t")
        if len(fields) != 2 or fields[0] in values:
            fail(f"malformed or duplicate key in {path.name}")
        values[fields[0]] = fields[1]
    return values


def parse_manifest(path: Path) -> tuple[dict[str, str], dict[str, str]]:
    lines = path.read_text(encoding="ascii").splitlines()
    if len(lines) != 4 or lines[0] != "A255885-MANIFEST\t1":
        fail("manifest record count or magic mismatch")
    names = lines[1].split("\t")
    initial = lines[2].split("\t")
    final = lines[3].split("\t")
    if len(names) != 25 or len(initial) != len(names) or len(final) != len(names):
        fail("manifest column mismatch")
    return dict(zip(names, initial)), dict(zip(names, final))


def parse_certificate(
    path: Path,
) -> tuple[dict[str, str], dict[int, int | None], dict[int, tuple[int, list[int]]]]:
    lines = path.read_text(encoding="ascii").splitlines()
    if not lines or lines[0] != "A255885-INDEPENDENT-CERT\t1":
        fail("certificate magic/version mismatch")
    metadata: dict[str, str] = {}
    minima: dict[int, int | None] = {}
    witnesses: dict[int, tuple[int, list[int]]] = {}
    for line in lines[1:]:
        fields = line.split("\t")
        if fields[0] == "M":
            if len(fields) != 3:
                fail("malformed certificate minimum row")
            n = int(fields[1])
            if n in minima:
                fail("duplicate certificate minimum")
            minima[n] = None if fields[2] == "NA" else int(fields[2])
        elif fields[0] == "W":
            if len(fields) != 5:
                fail("malformed certificate witness row")
            n, base, count = map(int, fields[1:4])
            if n in witnesses:
                fail("duplicate certificate witness row")
            values = [int(value) for value in fields[4].split(",")]
            if count != len(values):
                fail(f"witness count field mismatch at n={n}")
            witnesses[n] = (base, values)
        else:
            if len(fields) != 2 or fields[0] in metadata or minima:
                fail("malformed or duplicate certificate metadata")
            metadata[fields[0]] = fields[1]
    return metadata, minima, witnesses


def prime_table(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"
    for p in range(2, math.isqrt(limit) + 1):
        if sieve[p]:
            start = p * p
            sieve[start : limit + 1 : p] = b"\x00" * (((limit - start) // p) + 1)
    return [p for p in range(2, limit + 1) if sieve[p]]


def factor(n: int, primes: list[int]) -> list[tuple[int, int]]:
    original = n
    output: list[tuple[int, int]] = []
    for p in primes:
        if p * p > n:
            break
        if n % p:
            continue
        exponent = 0
        while n % p == 0:
            n //= p
            exponent += 1
        output.append((p, exponent))
    if n > 1:
        output.append((n, 1))
    if math.prod(p**e for p, e in output) != original:
        fail(f"factorization reconstruction failed for {original}")
    return output


def validate_witnesses(
    expected: set[int], witnesses: dict[int, tuple[int, list[int]]]
) -> int:
    if set(witnesses) != expected:
        fail("certificate witness-target set mismatch")
    maximum = max(c for _, values in witnesses.values() for c in values)
    primes = prime_table(math.isqrt(maximum) + 1)
    total = 0
    for n in sorted(witnesses):
        base, values = witnesses[n]
        if len(values) != n or values != sorted(set(values)):
            fail(f"witness order/cardinality mismatch at n={n}")
        total += len(values)
        for c in values:
            if not 4 <= c < base:
                fail(f"witness range mismatch at n={n}, c={c}")
            factors = factor(c, primes)
            if len(factors) == 1 and factors[0] == (c, 1):
                fail(f"prime listed as witness at n={n}, c={c}")
            modulus = c * c
            if pow(base, c - 1, modulus) != 1:
                fail(f"defining congruence failed at n={n}, c={c}")
            phi = c
            for p, _ in factors:
                phi -= phi // p
            if pow(base, phi, modulus) != 1:
                fail(f"Euler-Wieferich check failed at n={n}, c={c}")
            for p, exponent in factors:
                local = p ** (2 * exponent)
                if p == 2:
                    if base % local != 1:
                        fail(f"2-adic check failed at n={n}, c={c}")
                else:
                    h = math.gcd(c - 1, p - 1)
                    if pow(base, h, p) != 1 or pow(base, p - 1, local) != 1:
                        fail(f"odd-prime cascade failed at n={n}, c={c}, p={p}")
    return total


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--certificate", type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parent.parent
    certificate = args.certificate.resolve() if args.certificate else root / "validation/certificate.tsv"

    header, results = parse_results(root / "results/results.tsv")
    required_header = {
        "state": "complete",
        "B": str(BOUND),
        "modulus_frontier": "1999999999",
        "base_frontier": str(BOUND),
        "index_cap": str(INDEX_CAP),
        "counter_saturation": "255",
    }
    if header != required_header:
        fail("results metadata mismatch")
    if sum(value is not None for value in results) != EXPECTED_NON_NA:
        fail("results populated-slot count mismatch")
    if next(i for i, value in enumerate(results, 1) if value is None) != 86:
        fail("first unresolved result index is not 86")
    if parse_public_minima(root / "results/minima.tsv") != results:
        fail("normalized minima/raw result mismatch")

    contiguous = parse_pairs(root / "data/b255885.txt")
    if [n for n, _ in contiguous] != list(range(1, CONTIGUOUS_END + 1)):
        fail("b-file continuity mismatch")
    if [value for _, value in contiguous] != results[:CONTIGUOUS_END]:
        fail("b-file/result identity mismatch")

    sparse = parse_pairs(root / "data/sparse.txt")
    if len(sparse) != EXPECTED_SPARSE or [n for n, _ in sparse] != sorted({n for n, _ in sparse}):
        fail("sparse table order/cardinality mismatch")
    expected_sparse = [
        (n, value)
        for n, value in enumerate(results, 1)
        if n > CONTIGUOUS_END and value is not None
    ]
    if sparse != expected_sparse:
        fail("sparse table/result identity mismatch")

    config = parse_key_values(root / "validation/config.tsv", "A255885-CONFIG\t1")
    if config.get("B") != str(BOUND) or config.get("segment_span") != "1999999996":
        fail("official configuration boundary mismatch")
    if config.get("source_sha256") != "eab592a8c6a3caf75fa868724765800347516b23b65cc843c6d0a82c205aa1a5":
        fail("official source identity mismatch")
    if config.get("makefile_sha256") != "5449bb5683720e9992cb89b73da410b2f9c908beea2c0d35d888393ceff0b095":
        fail("official Makefile identity mismatch")

    initial, final = parse_manifest(root / "validation/manifest.tsv")
    if initial["state"] != "incomplete" or initial["frontier"] != "3":
        fail("initial manifest state mismatch")
    required_final = {
        "segment": "1",
        "c_lo": "4",
        "c_hi": "1999999999",
        "frontier": "1999999999",
        "roots_generated": "273026083586",
        "increments": "628611669",
        "saturations": "0",
        "composites": "1901777711",
        "threads": "5",
        "state": "complete",
        "active_slot": "1",
    }
    if any(final.get(key) != value for key, value in required_final.items()):
        fail("final manifest content mismatch")

    metadata, minima, witnesses = parse_certificate(certificate)
    required_certificate = {
        "B": str(BOUND),
        "threads": "16",
        "roots": final["roots_generated"],
        "increments": final["increments"],
        "composites": final["composites"],
        "full_counter_comparison": "PASS",
        "minima_comparison": "PASS",
        "new_target_count": str(EXPECTED_NEW),
    }
    if any(metadata.get(key) != value for key, value in required_certificate.items()):
        fail("certificate metadata mismatch")
    expected_minima = {n: value for n, value in enumerate(results, 1)}
    if minima != expected_minima:
        fail("certificate minima/result mismatch")
    target_set = {n for n, value in expected_minima.items() if n >= 25 and value is not None}
    if len(target_set) != EXPECTED_NEW:
        fail("new-target count mismatch")
    for n, (base, _) in witnesses.items():
        if base != expected_minima[n]:
            fail(f"witness base/minimum mismatch at n={n}")
    witness_total = validate_witnesses(target_set, witnesses)
    if witness_total != EXPECTED_WITNESSES:
        fail("total witness count mismatch")

    print("DATA VALIDATION PASS")
    print(f"results_sha256={sha256(root / 'results/results.tsv')}")
    print(f"contiguous_rows={len(contiguous)} sparse_rows={len(sparse)}")
    print(f"populated_slots={EXPECTED_NON_NA} empty_slots={INDEX_CAP - EXPECTED_NON_NA}")
    print(f"first_unresolved=86 witness_checks={witness_total}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (AssertionError, OSError, ValueError) as error:
        print(f"DATA VALIDATION FAIL: {error}")
        raise SystemExit(1)
