#!/usr/bin/env bash
set -euo pipefail
PORT=${1:-8080}
DIR=$(cd "$(dirname "$0")" && pwd)
cd "$DIR"
python3 -m http.server "$PORT"