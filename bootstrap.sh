#!/usr/bin/env bash
set -euo pipefail
PROJECT_ROOT="${PROJECT_ROOT:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)}"
export PROJECT_ROOT
if ! command -v python3 >/dev/null 2>&1; then
  echo 'Bootstrap requires Python 3.10+ (standard library only), Bash, and writable project storage.' >&2
  echo 'Install Python using your OS/WSL or an approved HPC module, then rerun ./bootstrap.sh.' >&2
  exit 2
fi
exec python3 "$PROJECT_ROOT/scripts/bootstrap.py" "$@"
