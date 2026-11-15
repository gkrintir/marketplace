#!/usr/bin/env bash
# Install cms-ml/marketplace skills and custom subagents into Google Antigravity.
#
# Usage:
#   scripts/install_antigravity.sh [--symlink|--copy] [--dest ~/.gemini/config] [--dry-run]
#
# Stdout: KEY=VALUE summary (SKILLS_INSTALLED, AGENTS_INSTALLED, DEST, MODE, DRY_RUN).
# Stderr: progress diagnostics.
# Exit codes: 0 success, 2 invalid option or missing plugins/ directory.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
DEST_DIR="${HOME}/.gemini/config"
MODE="symlink"
DRY_RUN=0

usage() {
  cat >&2 <<'EOF'
Usage: install_antigravity.sh [OPTIONS]

Install CMS Marketplace Agent Skills and adapted custom subagents
(cms-analysis-coder, cms-results-validator) into Google Antigravity.

Options:
  --symlink         Symlink skill directories so git updates take effect immediately (default)
  --copy            Copy skill directories instead of symlinking
  --repo <dir>      Marketplace checkout root (default: parent of scripts/)
  --dest <dir>      Target Antigravity config root (default: $HOME/.gemini/config)
  --dry-run         Report planned actions without writing to disk
  -h, --help        Show this help message and exit

Exit codes:
  0  Success
  2  Invalid option or missing plugins/ directory
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --symlink)
      MODE="symlink"
      shift
      ;;
    --copy)
      MODE="copy"
      shift
      ;;
    --repo)
      [[ $# -ge 2 ]] || { echo "error: --repo requires a path" >&2; exit 2; }
      REPO_DIR="$2"
      shift 2
      ;;
    --dest)
      [[ $# -ge 2 ]] || { echo "error: --dest requires a path" >&2; exit 2; }
      DEST_DIR="$2"
      shift 2
      ;;
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "error: unknown option: $1" >&2
      usage
      exit 2
      ;;
  esac
done

if [[ ! -d "${REPO_DIR}/plugins" ]]; then
  echo "error: plugins/ directory not found under ${REPO_DIR}" >&2
  exit 2
fi

REPO_DIR="$(cd "${REPO_DIR}" && pwd)"
SKILLS_DEST="${DEST_DIR}/skills"
AGENTS_DEST="${DEST_DIR}/agents"

if [[ "${DRY_RUN}" -eq 0 ]]; then
  mkdir -p "${SKILLS_DEST}" "${AGENTS_DEST}"
fi

echo "Installing skills from ${REPO_DIR}/plugins into ${SKILLS_DEST} (${MODE}, dry_run=${DRY_RUN}) ..." >&2
skill_count=0
for skill_dir in "${REPO_DIR}"/plugins/*/skills/*; do
  [[ -f "${skill_dir}/SKILL.md" ]] || continue
  skill_name="$(basename "${skill_dir}")"
  if [[ "${DRY_RUN}" -eq 0 ]]; then
    rm -rf "${SKILLS_DEST}/${skill_name}"
    if [[ "${MODE}" == "symlink" ]]; then
      ln -sfn "${skill_dir}" "${SKILLS_DEST}/${skill_name}"
    else
      cp -R "${skill_dir}" "${SKILLS_DEST}/${skill_name}"
    fi
  fi
  skill_count=$((skill_count + 1))
done

echo "Installing adapted custom subagents into ${AGENTS_DEST} (dry_run=${DRY_RUN}) ..." >&2
agent_count="$(
  python3 - "${REPO_DIR}" "${AGENTS_DEST}" "${DRY_RUN}" <<'PY'
import glob
import os
import sys

repo_dir, agents_dest, dry_run = sys.argv[1], sys.argv[2], int(sys.argv[3])
count = 0
for src in sorted(glob.glob(os.path.join(repo_dir, "plugins", "*", "agents", "*.md"))):
    with open(src, encoding="utf-8") as f:
        text = f.read()
    if not text.startswith("---\n"):
        continue
    _, fm, body = text.split("---\n", 2)
    kept_lines = []
    skip_continuation = False
    for line in fm.splitlines():
        if line and not line[0].isspace():
            key = line.split(":", 1)[0].strip()
            skip_continuation = key not in ("name", "description")
        if not skip_continuation:
            kept_lines.append(line)
    kept_lines.extend([
        "mainAgent: true",
        "subagent: true",
        "model: inherit",
        "enable_write_tools: true",
        "tools:",
        "  - send_message",
        "  - view_file",
        "  - google_search",
        "  - replace_file_content",
        "  - write_to_file",
        "  - run_command",
        "  - manage_task",
        "  - list_skills",
        "  - lookup_skill",
        "  - lookup_skills",
        "  - get_skill_resources",
    ])
    if not dry_run:
        out_path = os.path.join(agents_dest, os.path.basename(src))
        with open(out_path, "w", encoding="utf-8") as out:
            out.write("---\n" + "\n".join(kept_lines) + "\n---\n" + body)
    count += 1
print(count)
PY
)"

echo "Done. Start a new Antigravity conversation (or run '/skills reload') to load installed items." >&2
echo "SKILLS_INSTALLED=${skill_count}"
echo "AGENTS_INSTALLED=${agent_count}"
echo "DEST=${DEST_DIR}"
echo "MODE=${MODE}"
echo "DRY_RUN=${DRY_RUN}"
