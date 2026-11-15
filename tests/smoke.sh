#!/usr/bin/env bash
# Syntax and lint check only. This does not test the analysis itself.
set -euo pipefail
cd "$(dirname "$0")/../scripts"
fail=0
for f in *.sh; do bash -n "$f" || { echo "SYNTAX FAIL $f"; fail=1; }; done
for f in *.py; do [ -e "$f" ] && { python3 -I -m py_compile "$f" || fail=1; }; done
command -v shellcheck >/dev/null && shellcheck -S warning *.sh || true
exit $fail
