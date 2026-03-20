#!/bin/bash
set -euo pipefail

if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <user@host> <remote_dir>"
  exit 1
fi

PROJECT_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
REMOTE_HOST="$1"
REMOTE_DIR="$2"

rsync -az --delete \
  --exclude .git \
  --exclude build \
  --exclude build-uppmax \
  --exclude figures \
  --exclude data/outputs \
  "$PROJECT_ROOT/" "$REMOTE_HOST:$REMOTE_DIR/"
