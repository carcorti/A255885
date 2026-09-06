#!/usr/bin/env python3
"""Independent, fail-closed small-range validator for the A255885 C program.

This validator deliberately uses direct Python modular exponentiation rather
than the production primitive-root/CRT kernel.  It runs only ephemeral,
explicitly non-campaign computations with B <= 2024, while also checking the
v6 result and saturation boundaries.
"""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import itertools
import os
from pathlib import Path
import signal
import shutil
import struct
import subprocess
import tempfile
from typing import Iterable

KNOWN = [17, 65, 145, 485, 649, 1297, 577, 2024]
INDEX_CAP = 254
SAT = 255
MAX_B = 2_000_000_000


def fail(message: str) -> None:
    raise AssertionError(message)


def run(
    args: list[str],
    *,
    env: dict[str, str] | None = None,
    expect_ok: bool = True,
) -> subprocess.CompletedProcess[str]:
    proc = subprocess.run(args, text=True, capture_output=True, env=env, check=False)
    if expect_ok and proc.returncode != 0:
        fail(f"command failed: {args}\nstdout={proc.stdout}\nstderr={proc.stderr}")
    if not expect_ok and proc.returncode == 0:
        fail(f"command unexpectedly succeeded: {args}")
    return proc


def sieve_composite(bound: int) -> list[bool]:
    prime = bytearray(b"\x01") * (bound + 1)
    if bound >= 0:
        prime[0] = 0
    if bound >= 1:
        prime[1] = 0
    for p in range(2, int(bound**0.5) + 1):
        if prime[p]:
            prime[p * p : bound + 1 : p] = b"\x00" * (((bound - p * p) // p) + 1)
    return [n >= 4 and not bool(prime[n]) for n in range(bound + 1)]


def direct_counts(bound: int) -> tuple[bytearray, list[list[int]]]:
    composite = sieve_composite(bound)
    counts = bytearray(bound + 1)
    witnesses: list[list[int]] = [[] for _ in range(bound + 1)]
    for c in range(4, bound):
        if not composite[c]:
            continue
        c2 = c * c
        for b in range(c + 1, bound + 1):
            if pow(b, c - 1, c2) == 1:
                witnesses[b].append(c)
                counts[b] = min(SAT, counts[b] + 1)
    return counts, witnesses


def factor(n: int) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            e = 0
            while n % p == 0:
                n //= p
                e += 1
            out.append((p, e))
        p += 1 if p == 2 else 2
    if n > 1:
        out.append((n, 1))
    return out


def phi(n: int) -> int:
    result = n
    for p, _ in factor(n):
        result -= result // p
    return result


def valuation(n: int, p: int) -> int:
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def kernel_formula(c: int) -> int:
    result = 1
    for p, _ in factor(c):
        if p != 2:
            result *= __import__("math").gcd(c - 1, p - 1)
    return result


def root_and_cascade_regressions() -> None:
    composite = sieve_composite(100)
    for c in range(4, 81):
        if not composite[c]:
            continue
        roots = {r for r in range(c * c) if pow(r, c - 1, c * c) == 1}
        if len(roots) != kernel_formula(c):
            fail(f"kernel cardinality mismatch at c={c}")
        local_sets = []
        for p, e in factor(c):
            m = p ** (2 * e)
            local_sets.append((m, {r for r in range(m) if pow(r, c - 1, m) == 1}))
        reconstructed = {
            r
            for r in range(c * c)
            if all((r % m) in local for m, local in local_sets)
        }
        if reconstructed != roots:
            fail(f"CRT root-set mismatch at c={c}")

    roots15 = {r for r in range(225) if pow(r, 14, 225) == 1}
    if len(roots15) != 4:
        fail("mandatory c=15 regression did not produce four roots")
    if __import__("math").gcd(14, phi(15)) != 2:
        fail("c=15 false-formula counterexample was not reproduced")

    for c in range(4, 101):
        if not composite[c]:
            continue
        k = c - 1
        for p, e in factor(c):
            m = p ** (2 * e)
            for b in range(1, 201):
                direct = pow(b, k, m) == 1
                if p == 2:
                    cascade = b % m == 1
                else:
                    h = __import__("math").gcd(k, p - 1)
                    cascade = pow(b, h, p) == 1 and pow(b, p - 1, m) == 1
                if direct != cascade:
                    fail(f"cascade mismatch at c={c}, p={p}, b={b}")


def directory_digest(path: Path) -> str:
    h = hashlib.sha256()
    for item in sorted(path.iterdir(), key=lambda p: p.name):
        h.update(item.name.encode())
        h.update(item.read_bytes())
    return h.hexdigest()


CHECKSUM_FILES = [
    "config.tsv",
    "cnt0.bin",
    "cnt1.bin",
    "manifest.tsv",
    "results.tsv",
    "spf.bin",
]


def rewrite_checksums(run_dir: Path) -> None:
    text = "".join(
        f"{hashlib.sha256((run_dir / name).read_bytes()).hexdigest()}  {name}\n"
        for name in CHECKSUM_FILES
    )
    (run_dir / "checksums.sha256").write_text(text, encoding="ascii")


def rewrite_counter(run_dir: Path, slot: int, payload: bytearray) -> str:
    path = run_dir / f"cnt{slot}.bin"
    data = bytearray(path.read_bytes())
    if len(data) != 64 + len(payload):
        fail("synthetic checkpoint size mismatch")
    data[40:64] = hashlib.sha256(payload).digest()[:24]
    data[64:] = payload
    path.write_bytes(data)
    return hashlib.sha256(data).hexdigest()


def mutate_last_manifest(run_dir: Path, changes: dict[int, str]) -> None:
    path = run_dir / "manifest.tsv"
    lines = path.read_text(encoding="ascii").splitlines()
    fields = lines[-1].split("\t")
    if len(fields) != 25:
        fail("manifest mutation fixture is malformed")
    for index, value in changes.items():
        fields[index] = value
    lines[-1] = "\t".join(fields)
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def expect_next_immutable_failure(binary: Path, run_dir: Path) -> None:
    before = directory_digest(run_dir)
    run([str(binary), "next", str(run_dir)], expect_ok=False)
    if directory_digest(run_dir) != before:
        fail("rejected persistence corruption modified live artifacts")


def read_run_counts(run_dir: Path) -> bytearray:
    lines = (run_dir / "manifest.tsv").read_text(encoding="ascii").splitlines()
    fields = lines[-1].split("\t")
    if len(fields) != 25:
        fail("validator saw malformed final manifest")
    slot = int(fields[18])
    B = int((run_dir / "config.tsv").read_text(encoding="ascii").split("B\t", 1)[1].splitlines()[0])
    data = (run_dir / f"cnt{slot}.bin").read_bytes()
    if len(data) != 64 + B + 1:
        fail("counter checkpoint size mismatch")
    if data[: len(b"A255885-BINARY")] != b"A255885-BINARY":
        fail("counter checkpoint magic mismatch")
    version, kind, header_B, elem = struct.unpack_from("<IIII", data, 16)
    count = struct.unpack_from("<Q", data, 32)[0]
    if (version, kind, header_B, elem, count) != (1, 2, B, 1, B + 1):
        fail("counter checkpoint header mismatch")
    payload = bytearray(data[64:])
    if any(v > SAT for v in payload):
        fail("counter saturation domain mismatch")
    return payload


def complete_run(
    binary: Path,
    root: Path,
    name: str,
    bound: int,
    threads: int,
    span: int | None,
    pause_after: int = 0,
) -> tuple[Path, bytearray]:
    run_dir = root / name
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    if span is not None:
        env["A255885_NONCAMPAIGN"] = "1"
        env["A255885_TEST_SEGMENT_SPAN"] = str(span)
    run([str(binary), "init", str(run_dir), str(bound)], env=env)

    before = directory_digest(run_dir)
    run([str(binary), "status", str(run_dir)], env=env)
    after = directory_digest(run_dir)
    if before != after:
        fail("status modified the run directory")

    segments = 0
    while True:
        status = run([str(binary), "status", str(run_dir)], env=env).stdout
        state = dict(line.split("\t", 1) for line in status.splitlines())["state"]
        if state == "complete":
            break
        run([str(binary), "next", str(run_dir)], env=env)
        segments += 1
        if pause_after and segments == pause_after:
            run([str(binary), "verify", str(run_dir)], env=env)
            run([str(binary), "status", str(run_dir)], env=env)
    run([str(binary), "verify", str(run_dir)], env=env)
    return run_dir, read_run_counts(run_dir)


def theorem_55_certificate(b: int, c: int) -> None:
    if pow(b, phi(c), c * c) != 1:
        fail(f"Euler-Wieferich necessary condition failed for b={b}, c={c}")
    factors = factor(c)
    product = 1
    for p, _ in factors:
        product *= p - 1
    for p, e in factors:
        if p == 2:
            if b % 4 == 1:
                lam = valuation(b - 1, 2) - 1
            else:
                lam = valuation(b + 1, 2) - 1
        else:
            lam = valuation(pow(b, p - 1) - 1, p) - 1
        if valuation(product, p) + lam < e:
            fail(f"Theorem 5.5 inequality failed for b={b}, c={c}, p={p}")


def expect_corruption_failure(binary: Path, source: Path, root: Path, name: str, mutate) -> None:
    target = root / name
    shutil.copytree(source, target)
    mutate(target)
    run([str(binary), "status", str(target)], expect_ok=False)


def fault_injection(binary: Path, completed: Path, root: Path) -> None:
    run([str(binary), "init", str(completed), "100"], expect_ok=False)

    def corrupt_active(d: Path) -> None:
        fields = (d / "manifest.tsv").read_text(encoding="ascii").splitlines()[-1].split("\t")
        p = d / f"cnt{fields[18]}.bin"
        data = bytearray(p.read_bytes())
        data[-1] ^= 1
        p.write_bytes(data)

    expect_corruption_failure(binary, completed, root, "bad-count", corrupt_active)
    expect_corruption_failure(
        binary,
        completed,
        root,
        "bad-manifest",
        lambda d: (d / "manifest.tsv").write_bytes((d / "manifest.tsv").read_bytes()[:-1]),
    )
    expect_corruption_failure(
        binary,
        completed,
        root,
        "bad-checksum",
        lambda d: (d / "checksums.sha256").write_text("0" * 64 + "  config.tsv\n", encoding="ascii"),
    )
    expect_corruption_failure(
        binary,
        completed,
        root,
        "leftover-tmp",
        lambda d: (d / "cnt0.bin.tmp").write_bytes(b"partial"),
    )


def commit_fault_recovery(binary: Path, root: Path) -> None:
    expected, _ = direct_counts(300)
    fault_codes = {
        "before_rename": 91,
        "after_rename": 92,
        "after_manifest": 93,
        "after_results": 94,
    }
    for stage, expected_code in fault_codes.items():
        run_dir = root / f"fault-{stage}"
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = "2"
        env["A255885_NONCAMPAIGN"] = "1"
        env["A255885_TEST_SEGMENT_SPAN"] = "29"
        run([str(binary), "init", str(run_dir), "300"], env=env)
        fault_env = env.copy()
        fault_env["A255885_TEST_FAULT"] = stage
        fault = run([str(binary), "next", str(run_dir)], env=fault_env, expect_ok=False)
        if fault.returncode != expected_code or f"test_fault\tstage={stage}" not in fault.stderr:
            fail(f"fault stage {stage} did not exit at its exact boundary")
        run([str(binary), "status", str(run_dir)], env=env, expect_ok=False)
        recovered = run([str(binary), "next", str(run_dir)], env=env)
        if stage == "before_rename":
            recovery_marker = "recovery\tdiscarded=cnt"
        elif stage == "after_rename":
            recovery_marker = "recovery\tstate=discarded_uncommitted_checkpoint"
        else:
            recovery_marker = "recovery\tstate=completed_manifest_commit"
        if recovery_marker not in recovered.stderr:
            fail(f"fault stage {stage} used the wrong recovery lineage")
        while True:
            status = run([str(binary), "status", str(run_dir)], env=env).stdout
            state = dict(line.split("\t", 1) for line in status.splitlines())["state"]
            if state == "complete":
                break
            run([str(binary), "next", str(run_dir)], env=env)
        run([str(binary), "verify", str(run_dir)], env=env)
        if read_run_counts(run_dir) != expected:
            fail(f"commit-fault recovery changed scientific counts at {stage}")


def manifest_authority_regressions(binary: Path, root: Path) -> None:
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "2"
    env["A255885_NONCAMPAIGN"] = "1"
    env["A255885_TEST_SEGMENT_SPAN"] = "29"
    base = root / "authority-base"
    run([str(binary), "init", str(base), "300"], env=env)
    run([str(binary), "next", str(base)], env=env)
    one_row_manifest = (base / "manifest.tsv").read_bytes()

    cases: list[tuple[str, dict[int, str]]] = [
        ("forged-frontier", {2: "299", 4: "299", 17: "complete", 24: "complete"}),
        ("wrong-slot", {18: "0"}),
        ("wrong-state", {17: "complete"}),
        ("wrong-stop", {24: "complete"}),
    ]
    for name, changes in cases:
        target = root / name
        shutil.copytree(base, target)
        mutate_last_manifest(target, changes)
        expect_next_immutable_failure(binary, target)

    run([str(binary), "next", str(base)], env=env)
    stale = root / "stale-manifest"
    shutil.copytree(base, stale)
    (stale / "manifest.tsv").write_bytes(one_row_manifest)
    expect_next_immutable_failure(binary, stale)

    partial = root / "partial-manifest-temp"
    shutil.copytree(base, partial)
    (partial / "manifest.tsv.tmp").write_bytes(b"partial appended row")
    run([str(binary), "status", str(partial)], expect_ok=False)
    run([str(binary), "next", str(partial)], env=env)
    run([str(binary), "verify", str(partial)], env=env)

    bad_results = root / "bad-results-current-manifest"
    shutil.copytree(base, bad_results)
    (bad_results / "results.tsv").write_bytes(b"corrupt derived results\n")
    expect_next_immutable_failure(binary, bad_results)


def run_directory_exclusion_regression(binary: Path, root: Path) -> None:
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "2"
    env["A255885_NONCAMPAIGN"] = "1"
    env["A255885_TEST_SEGMENT_SPAN"] = "2020"

    rejected_hook = root / "rejected-init-lock-hook"
    rejected_env = os.environ.copy()
    rejected_env["A255885_TEST_STOP_AFTER_INIT_LOCK"] = "1"
    rejected = run(
        [str(binary), "init", str(rejected_hook), "2024"],
        env=rejected_env,
        expect_ok=False,
    )
    if "init-lock stop is restricted to recorded non-campaign tests" not in rejected.stderr:
        fail("init-lock test hook was not rejected outside the non-campaign gate")
    if rejected_hook.exists():
        fail("rejected init-lock test hook created a run directory")

    held = root / "held-run-directory-lock"
    run([str(binary), "init", str(held), "2024"], env=env)
    before = directory_digest(held)
    lock_fd = os.open(held, os.O_RDONLY | os.O_DIRECTORY)
    try:
        fcntl.flock(lock_fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        for command in ("status", "verify", "next"):
            blocked = run([str(binary), command, str(held)], env=env, expect_ok=False)
            if "run directory is already in use by another a255885 process" not in blocked.stderr:
                fail(f"held run-directory lock did not exclude {command} exactly")
        if directory_digest(held) != before:
            fail("lock rejection modified the run directory")
    finally:
        os.close(lock_fd)

    shared_fd = os.open(held, os.O_RDONLY | os.O_DIRECTORY)
    try:
        fcntl.flock(shared_fd, fcntl.LOCK_SH | fcntl.LOCK_NB)
        run([str(binary), "status", str(held)], env=env)
        run([str(binary), "verify", str(held)], env=env)
        before = directory_digest(held)
        blocked = run([str(binary), "next", str(held)], env=env, expect_ok=False)
        if "run directory is already in use by another a255885 process" not in blocked.stderr:
            fail("shared reader lock did not exclude next exactly")
        if directory_digest(held) != before:
            fail("writer rejection under a shared reader lock modified the run directory")
    finally:
        os.close(shared_fd)

    init_race = root / "init-next-exclusion"
    init_env = env.copy()
    init_env["A255885_TEST_STOP_AFTER_INIT_LOCK"] = "1"
    init_proc = subprocess.Popen(
        [str(binary), "init", str(init_race), "2024"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=init_env,
    )
    stopped = False
    try:
        waited_pid, wait_status = os.waitpid(init_proc.pid, os.WUNTRACED)
        if waited_pid != init_proc.pid or not os.WIFSTOPPED(wait_status):
            fail("init-lock fixture did not reach its deterministic stop")
        stopped = True
        before = directory_digest(init_race)
        blocked = run([str(binary), "next", str(init_race)], env=env, expect_ok=False)
        if "run directory is already in use by another a255885 process" not in blocked.stderr:
            fail("init did not exclude next while creating the persistent state")
        if directory_digest(init_race) != before:
            fail("next rejection modified the stopped init directory")
        os.kill(init_proc.pid, signal.SIGCONT)
        stopped = False
        stdout, stderr = init_proc.communicate()
        if init_proc.returncode != 0:
            fail(f"stopped init did not resume cleanly\nstdout={stdout}\nstderr={stderr}")
    finally:
        if stopped:
            os.kill(init_proc.pid, signal.SIGCONT)
        if init_proc.poll() is None:
            init_proc.kill()
            init_proc.wait()
    run([str(binary), "verify", str(init_race)], env=env)

    concurrent = root / "concurrent-next"
    run([str(binary), "init", str(concurrent), "2024"], env=env)
    wrappers: list[subprocess.Popen[str]] = []
    shell = 'IFS= read -r ready; exec "$1" next "$2"'
    for _ in range(16):
        wrappers.append(
            subprocess.Popen(
                ["/bin/sh", "-c", shell, "a255885-lock-test", str(binary), str(concurrent)],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                env=env,
            )
        )
    for proc in wrappers:
        if proc.stdin is None:
            fail("concurrent lock fixture has no stdin barrier")
        proc.stdin.write("go\n")
        proc.stdin.close()

    outcomes: list[tuple[int, str, str]] = []
    for proc in wrappers:
        stdout = proc.stdout.read() if proc.stdout is not None else ""
        stderr = proc.stderr.read() if proc.stderr is not None else ""
        outcomes.append((proc.wait(), stdout, stderr))
    successes = [item for item in outcomes if item[0] == 0]
    failures = [item for item in outcomes if item[0] != 0]
    if len(successes) != 1 or len(failures) != 15:
        fail(f"concurrent next exclusion count mismatch: success={len(successes)} failure={len(failures)}")
    if any("run directory is already in use by another a255885 process" not in item[2]
           for item in failures):
        fail("concurrent next used a failure path other than run-directory exclusion")
    manifest_lines = (concurrent / "manifest.tsv").read_text(encoding="ascii").splitlines()
    if len(manifest_lines) != 4 or manifest_lines[-1].split("\t")[0] != "1":
        fail("concurrent next committed other than exactly one segment")
    run([str(binary), "verify", str(concurrent)], env=env)
    direct, _ = direct_counts(2024)
    if read_run_counts(concurrent) != direct:
        fail("run-directory exclusion changed scientific counts")


def segment_and_thread_regressions(binary: Path, root: Path) -> None:
    expected, _ = direct_counts(300)
    explicit, counts = complete_run(binary, root, "explicit-span", 300, 2, None)
    if counts != expected:
        fail("default segment result changed")

    configured = root / "configured-span"
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "2"
    run([str(binary), "init", str(configured), "300", "29"], env=env)
    while True:
        status = run([str(binary), "status", str(configured)], env=env).stdout
        if dict(line.split("\t", 1) for line in status.splitlines())["state"] == "complete":
            break
        run([str(binary), "next", str(configured)], env=env)
    verify = run([str(binary), "verify", str(configured)], env=env).stdout
    if "scope=identity+persistence+coverage" not in verify or "arithmetic=not_recomputed" not in verify:
        fail("verify scope is not explicit")
    if read_run_counts(configured) != expected:
        fail("explicit segment span changed scientific counts")

    for name, span in (("span-zero", "0"), ("span-large", "297")):
        rejected = root / name
        run([str(binary), "init", str(rejected), "300", span], expect_ok=False)
        if rejected.exists():
            fail("invalid segment span left a run directory")

    thread_dir = root / "actual-threads"
    thread_env = os.environ.copy()
    thread_env["OMP_NUM_THREADS"] = "8"
    thread_env["OMP_THREAD_LIMIT"] = "2"
    run([str(binary), "init", str(thread_dir), "100", "96"], env=thread_env)
    run([str(binary), "next", str(thread_dir)], env=thread_env)
    fields = (thread_dir / "manifest.tsv").read_text(encoding="ascii").splitlines()[-1].split("\t")
    if fields[16] != "2":
        fail(f"manifest recorded requested rather than actual threads: {fields[16]}")


def saturation_telemetry_regression(binary: Path, root: Path) -> None:
    for initial, expected, expected_saturations in (
        (253, 254, "0"),
        (254, 255, "1"),
        (255, 255, "0"),
    ):
        seed = root / f"saturation-seed-{initial}"
        env = os.environ.copy()
        env["A255885_NONCAMPAIGN"] = "1"
        env["A255885_TEST_SEGMENT_SPAN"] = "13"
        run([str(binary), "init", str(seed), "17"], env=env)
        payload = read_run_counts(seed)
        payload[17] = initial
        count_hash = rewrite_counter(seed, 0, payload)
        mutate_last_manifest(seed, {21: count_hash})
        rewrite_checksums(seed)

        outputs: list[tuple[bytearray, str]] = []
        for threads in (1, 4):
            target = root / f"saturation-{initial}-t{threads}"
            shutil.copytree(seed, target)
            thread_env = os.environ.copy()
            thread_env["OMP_NUM_THREADS"] = str(threads)
            run([str(binary), "next", str(target)], env=thread_env)
            run([str(binary), "verify", str(target)], env=thread_env)
            fields = (target / "manifest.tsv").read_text(encoding="ascii").splitlines()[-1].split("\t")
            outputs.append((read_run_counts(target), fields[13]))
        if outputs[0][0] != outputs[1][0]:
            fail(f"saturation payload is thread-dependent from initial={initial}")
        if any(item[0][17] != expected or item[1] != expected_saturations for item in outputs):
            fail(f"saturation boundary mismatch from initial={initial}")


def cli_matrix(binary: Path, root: Path) -> None:
    bad = [
        [str(binary)],
        [str(binary), "unknown"],
        [str(binary), "init", str(root / "x")],
        [str(binary), "init", str(root / "zero"), "0"],
        [str(binary), "init", str(root / "neg"), "-1"],
        [str(binary), "init", str(root / "plus"), "+2"],
        [str(binary), "init", str(root / "space"), " 2"],
        [str(binary), "init", str(root / "over"), str(MAX_B + 1)],
        [str(binary), "status"],
        [str(binary), "next", str(root / "missing"), "extra"],
    ]
    for args in bad:
        run(args, expect_ok=False)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    args = parser.parse_args()
    binary = args.binary.resolve()
    if not binary.is_file():
        fail("binary does not exist")

    root_and_cascade_regressions()
    direct, witnesses = direct_counts(2024)
    minima = []
    for target in range(1, 9):
        found = next((b for b in range(2, 2025) if direct[b] == target), None)
        minima.append(found)
    if minima != KNOWN:
        fail(f"direct replay mismatch: {minima}")

    for target, b in enumerate(KNOWN, start=1):
        if len(witnesses[b]) != target:
            fail(f"witness cardinality mismatch at a({target})")
        for c in witnesses[b]:
            if not (2 <= c < b) or pow(b, c - 1, c * c) != 1:
                fail("direct witness reconstruction failed")
            theorem_55_certificate(b, c)

    with tempfile.TemporaryDirectory(prefix="a255885-validation-") as tmp:
        root = Path(tmp)
        run1, c1 = complete_run(binary, root, "t1", 2024, 1, None)
        _, c2 = complete_run(binary, root, "t2", 2024, 2, 37, pause_after=5)
        _, c4 = complete_run(binary, root, "t4", 2024, 4, 251)
        if c1 != direct or c2 != direct or c4 != direct:
            fail("C/Python, thread-count, or segmentation equality failed")

        results = (run1 / "results.tsv").read_text(encoding="ascii").splitlines()
        header = dict(line.split("\t", 1) for line in results[1:7])
        if header.get("index_cap") != str(INDEX_CAP) or header.get("counter_saturation") != str(SAT):
            fail("result-domain header mismatch")
        rows = results[8:]
        if len(rows) != INDEX_CAP:
            fail(f"result row count mismatch: {len(rows)}")
        reported: list[int | None] = []
        for target, line in enumerate(rows, start=1):
            index_text, base_text = line.split("\t")
            if int(index_text) != target:
                fail("result index order mismatch")
            reported.append(None if base_text == "NA" else int(base_text))
        expected_minima = [
            next((b for b in range(2, 2025) if direct[b] == target), None)
            for target in range(1, INDEX_CAP + 1)
        ]
        if reported != expected_minima:
            fail("single-pass result minima differ from independent Python minima")
        if reported[:8] != KNOWN:
            fail(f"C known-term replay mismatch: {reported[:8]}")

        fault_injection(binary, run1, root)
        commit_fault_recovery(binary, root)
        manifest_authority_regressions(binary, root)
        run_directory_exclusion_regression(binary, root)
        segment_and_thread_regressions(binary, root)
        saturation_telemetry_regression(binary, root)
        cli_matrix(binary, root)

    print("VALIDATION PASS")
    print("direct_definition_B=2024: PASS")
    print("known_minima_a1_a8: 17,65,145,485,649,1297,577,2024")
    print("c15_exact_roots=4_false_gcd_formula=2: PASS")
    print("local_full_CRT_root_sets_c<=80: PASS")
    print("primewise_cascade_c<=100_b<=200: PASS")
    print("C_vs_independent_Python_counts_B=2024: PASS")
    print("threads_1_2_4_byte_identity: PASS")
    print("segment_spans_default_37_251_and_resume: PASS")
    print("witness_lists_and_Theorem_5_5_certificates: PASS")
    print("status_read_only_and_verify_fail_closed: PASS")
    print("fault_injection_and_CLI_matrix: PASS")
    print("commit_faults_before_after_rename_after_manifest_and_after_results: PASS")
    print("manifest_transition_forgery_stale_and_temp_recovery: PASS")
    print("run_directory_exclusion_writers_readers_init_and_preheld_lock: PASS")
    print("explicit_segment_span_and_actual_thread_telemetry: PASS")
    print("authoritative_saturation_transition_thread_independence: PASS")
    print("single_pass_minima_rows_1_254_vs_independent_Python: PASS")
    print("verify_scope_identity_persistence_coverage_only: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
