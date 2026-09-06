#!/bin/sh
set -eu

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
MODE=${1:-quick}

cleanup_builds() {
    make -C "$ROOT/src" clean >/dev/null 2>&1 || true
    make -C "$ROOT/validation" clean >/dev/null 2>&1 || true
}

quick() {
    trap cleanup_builds EXIT HUP INT TERM
    python3 "$ROOT/validation/check_package.py"
    python3 "$ROOT/validation/check_data.py"
    make -C "$ROOT/src" test
    make -C "$ROOT/validation" test
    cleanup_builds
    python3 "$ROOT/validation/check_package.py"
    echo "PUBLIC VALIDATION PASS"
}

paper() {
    if ! command -v pdflatex >/dev/null 2>&1; then
        echo "pdflatex is required for paper mode" >&2
        exit 1
    fi
    TMP=$(mktemp -d "${TMPDIR:-/tmp}/a255885-paper.XXXXXX")
    cleanup_paper() {
        rm -f -- "$TMP"/*
        rmdir -- "$TMP"
    }
    trap cleanup_paper EXIT HUP INT TERM
    pdflatex -draftmode -interaction=nonstopmode -halt-on-error \
        -output-directory="$TMP" "$ROOT/paper/A255885_v3.tex" >/dev/null
    pdflatex -draftmode -interaction=nonstopmode -halt-on-error \
        -output-directory="$TMP" "$ROOT/paper/A255885_v3.tex" >/dev/null
    if grep -E "Warning|undefined|Overfull|Underfull|Error" "$TMP/A255885_v3.log"; then
        echo "PAPER VALIDATION FAIL: log warning or error" >&2
        exit 1
    fi
    test ! -e "$ROOT/paper/A255885_v3.pdf"
    echo "PAPER DRAFTMODE PASS"
}

full() {
    if [ "$#" -ne 1 ] || [ "${A255885_ALLOW_FULL:-}" != "YES" ]; then
        echo "Usage: A255885_ALLOW_FULL=YES sh validation/run.sh full /absolute/new/run-dir" >&2
        exit 1
    fi
    RUN_DIR=$1
    case "$RUN_DIR" in
        /*) ;;
        *) echo "full mode requires an absolute run path" >&2; exit 1 ;;
    esac
    if [ -e "$RUN_DIR" ]; then
        echo "full mode refuses an existing run path: $RUN_DIR" >&2
        exit 1
    fi
    PARENT=$(dirname -- "$RUN_DIR")
    if [ ! -d "$PARENT" ]; then
        echo "run parent does not exist: $PARENT" >&2
        exit 1
    fi
    AVAILABLE_KIB=$(awk '/^MemAvailable:/ {print $2}' /proc/meminfo)
    if [ "${AVAILABLE_KIB:-0}" -lt 31457280 ]; then
        echo "full mode requires at least 31457280 KiB MemAvailable" >&2
        exit 1
    fi
    AVAILABLE_KB=$(df -Pk "$PARENT" | awk 'NR==2 {print $4}')
    if [ "${AVAILABLE_KB:-0}" -lt 16000000 ]; then
        echo "full mode requires at least 16 GB free disk space" >&2
        exit 1
    fi
    trap cleanup_builds EXIT HUP INT TERM
    make -C "$ROOT/src" release
    make -C "$ROOT/validation" release
    "$ROOT/src/a255885" init "$RUN_DIR" 2000000000 1999999996
    OMP_NUM_THREADS=5 OMP_DYNAMIC=FALSE OMP_THREAD_LIMIT=5 \
        "$ROOT/src/a255885" next "$RUN_DIR"
    "$ROOT/src/a255885" verify "$RUN_DIR"
    cmp "$ROOT/results/results.tsv" "$RUN_DIR/results.tsv"
    SLOT=$(awk -F '\t' 'END {print $19}' "$RUN_DIR/manifest.tsv")
    CERT="$RUN_DIR/replay_cert.tsv"
    OMP_NUM_THREADS=16 OMP_DYNAMIC=FALSE OMP_THREAD_LIMIT=16 \
        "$ROOT/validation/a255885-independent" validate \
        "$RUN_DIR/results.tsv" "$RUN_DIR/cnt${SLOT}.bin" "$CERT" --threads 16
    python3 "$ROOT/validation/check_data.py" --certificate "$CERT"
    cleanup_builds
    echo "FULL FRONTIER REPRODUCTION PASS"
    echo "retained_run_dir=$RUN_DIR"
}

case "$MODE" in
    quick)
        [ "$#" -eq 1 ] || { echo "Usage: sh validation/run.sh quick" >&2; exit 1; }
        quick
        ;;
    paper)
        [ "$#" -eq 1 ] || { echo "Usage: sh validation/run.sh paper" >&2; exit 1; }
        paper
        ;;
    full)
        shift
        full "$@"
        ;;
    *)
        echo "Usage: sh validation/run.sh {quick|paper|full}" >&2
        exit 1
        ;;
esac

