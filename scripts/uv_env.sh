#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
UV_BIN="${REPO_ROOT}/.tools/uv"
VENV_DIR="${REPO_ROOT}/.venv"
CACHE_DIR="${REPO_ROOT}/.cache"

# In some sandboxed environments $HOME may not be writable. Force uv/xdg caches into the repo.
export XDG_CACHE_HOME="${CACHE_DIR}"
export UV_CACHE_DIR="${CACHE_DIR}/uv"

usage() {
  cat <<EOF
Usage:
  scripts/uv_env.sh install-uv
  scripts/uv_env.sh venv
  scripts/uv_env.sh run-dml

Notes:
  - Installs uv into ${UV_BIN}
  - Creates a venv at ${VENV_DIR} using system site-packages (so numpy/pandas can be reused)
EOF
}

install_uv() {
  mkdir -p "${REPO_ROOT}/.tools"
  mkdir -p "${UV_CACHE_DIR}"
  if [[ -x "${UV_BIN}" ]]; then
    echo "uv already present at ${UV_BIN}"
    exit 0
  fi

  # Download the official installer script and install uv into .tools/
  # Requires network access.
  curl -LsSf https://astral.sh/uv/install.sh | UV_INSTALL_DIR="${REPO_ROOT}/.tools" sh
  test -x "${UV_BIN}"
  "${UV_BIN}" --version
}

make_venv() {
  if [[ ! -x "${UV_BIN}" ]]; then
    echo "uv not found at ${UV_BIN}. Run: scripts/uv_env.sh install-uv" >&2
    exit 1
  fi

  mkdir -p "${UV_CACHE_DIR}"
  "${UV_BIN}" venv --system-site-packages "${VENV_DIR}"
  echo "Created venv at ${VENV_DIR}"
  echo "Python: ${VENV_DIR}/bin/python"
}

run_dml() {
  if [[ ! -x "${VENV_DIR}/bin/python" ]]; then
    echo "venv not found at ${VENV_DIR}. Run: scripts/uv_env.sh venv" >&2
    exit 1
  fi
  "${VENV_DIR}/bin/python" "${REPO_ROOT}/python/02_dml_example.py"
}

cmd="${1:-}"
case "${cmd}" in
  install-uv) install_uv ;;
  venv) make_venv ;;
  run-dml) run_dml ;;
  *) usage; exit 2 ;;
esac
