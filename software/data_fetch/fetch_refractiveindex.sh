#!/usr/bin/env bash
# Клон/обновление базы refractiveindex.info (CC0) в data/reference/optical_constants/.
# Затем — индекс файлов для материалов из materials.tsv.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
need git jq

repo_url="https://github.com/polyanskiy/refractiveindex.info-database.git"
dest="$SSP_REF/optical_constants/refractiveindex.info-database"
mkdir -p "$(dirname "$dest")"

if [[ -d "$dest/.git" ]]; then
  log "обновление $dest"
  git -C "$dest" pull --ff-only --depth 1
else
  log "клонирование $repo_url"
  git clone --depth 1 "$repo_url" "$dest"
fi
commit="$(git -C "$dest" rev-parse HEAD)"

write_meta "$dest.meta.json" "refractiveindex.info database" "$repo_url" "CC0-1.0" \
  "M. N. Polyanskiy, Refractiveindex.info database of optical constants, Sci. Data 11, 94 (2024)" \
  "Каждый YAML-файл — набор конкретных авторов; ссылка на первоисточник внутри файла (REFERENCES)." \
  "" "git $commit"

# Индекс нужных материалов: material<TAB>относительный путь<TAB>первая строка REFERENCES
index="$SSP_REF/optical_constants/index_starter_set.tsv"
data_root="$dest/database/data"
printf 'material\tfile\treference\n' > "$index"
materials | while IFS=$'\t' read -r name _ _ glob; do
  shopt -s nullglob
  files=( "$data_root"/$glob )
  shopt -u nullglob
  [[ ${#files[@]} -eq 0 ]] && log "RII: не найдено для $name ($glob)"
  for f in "${files[@]}"; do
    ref="$(grep -m1 -A1 '^REFERENCES' "$f" | tail -n1 | sed 's/^[[:space:]"]*//; s/<[^>]*>//g' | cut -c1-200)"
    printf '%s\t%s\t%s\n' "$name" "${f#$dest/}" "$ref" >> "$index"
  done
done
log "индекс: $index ($(($(wc -l < "$index") - 1)) файлов)"
