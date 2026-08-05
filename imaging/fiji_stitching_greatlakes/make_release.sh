#!/usr/bin/env bash
set -euo pipefail
BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
DESTINATION="${1:-${BUNDLE_DIR}/dist}"
mkdir -p "${DESTINATION}"
ARCHIVE="${DESTINATION}/fiji_stitching_greatlakes_$(date +%Y%m%d).tar.gz"
tar --exclude='./dist' -C "${BUNDLE_DIR}" -czf "${ARCHIVE}" .
sha256sum "${ARCHIVE}" > "${ARCHIVE}.sha256"
echo "${ARCHIVE}"
echo "${ARCHIVE}.sha256"
