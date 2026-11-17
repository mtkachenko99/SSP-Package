# Общие функции для скриптов выгрузки (подключается через source).
set -euo pipefail

SSP_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SSP_REF="$SSP_ROOT/data/reference"
SSP_MATERIALS="${SSP_MATERIALS:-$SSP_ROOT/software/data_fetch/materials.tsv}"
SSP_UA="ssp-nano-data-fetch/0.1"

log() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }

need() {
  local c; for c in "$@"; do
    command -v "$c" >/dev/null || { log "нет программы: $c"; exit 1; }
  done
}

# Строки materials.tsv без комментариев: name cod_formula jarvis_formula rii_glob
materials() { grep -v '^#' "$SSP_MATERIALS" | grep -v '^[[:space:]]*$'; }

# write_meta <файл.meta.json> <source> <url> <license> <citation> [notes] [level_of_theory] [version]
# Метаданные происхождения по методологии M2.2.
write_meta() {
  jq -n --arg source "$2" --arg url "$3" --arg license "$4" --arg citation "$5" \
        --arg notes "${6:-}" --arg lot "${7:-}" --arg version "${8:-}" \
        --arg retrieved "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
        '{source:$source, url:$url, retrieved:$retrieved, version:$version, license:$license,
          citation:$citation, level_of_theory:$lot, notes:$notes}' > "$1"
}
