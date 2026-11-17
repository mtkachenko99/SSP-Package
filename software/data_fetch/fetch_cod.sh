#!/usr/bin/env bash
# Выгрузка CIF из Crystallography Open Database (CC0) по формулам из materials.tsv.
# Для каждого материала: index.json (результат поиска COD со всеми полями, включая R-факторы),
# <codid>.cif и index.tsv (краткая сводка для фильтрации).
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
need curl jq

base="https://www.crystallography.net/cod"
max_per_material="${COD_MAX:-300}"   # ограничение числа CIF на материал (SiO2 даёт сотни записей)
COD_DELAY="${COD_DELAY:-1.5}"        # пауза между запросами, с: COD сбрасывает соединения при частых запросах
CURL="curl -sSf --max-time 120 --retry 5 --retry-delay 10 --retry-all-errors -A $SSP_UA"

materials | while IFS=$'\t' read -r name formula _ _; do
  dir="$SSP_REF/structures/cod/$name"
  mkdir -p "$dir"
  q="$(jq -rn --arg f "$formula" '$f|@uri')"
  url="$base/result?formula=$q&format=json"
  log "COD: $name ($formula)"
  if ! $CURL "$url" -o "$dir/index.json.part"; then
    log "  поиск не удался (сервер недоступен?) — пропуск, перезапустите позже"; rm -f "$dir/index.json.part"; continue
  fi
  mv "$dir/index.json.part" "$dir/index.json"
  n="$(jq length "$dir/index.json")"

  jq -r '["codid","sg","sgNumber","a","b","c","alpha","beta","gamma","Robs","celltemp","cellpressure","mineral","chemname","year","doi"],
         (.[] | [.file,.sg,.sgNumber,.a,.b,.c,.alpha,.beta,.gamma,.Robs,.celltemp,.cellpressure,.mineral,.chemname,.year,.doi])
         | map(. // "") | @tsv' "$dir/index.json" > "$dir/index.tsv"

  # Сначала записи с наименьшим R-фактором (без R — в конце)
  jq -r 'sort_by(.Robs // "9" | tonumber? // 9) | .[].file' "$dir/index.json" \
    | head -n "$max_per_material" | while read -r id; do
      [[ -s "$dir/$id.cif" ]] && continue
      $CURL "$base/$id.cif" -o "$dir/$id.cif" \
        || { log "  не удалось: $id"; rm -f "$dir/$id.cif"; }
      sleep "$COD_DELAY"
    done
  got="$(find "$dir" -name '*.cif' | wc -l)"
  log "  записей: $n, CIF скачано: $got"

  write_meta "$dir.meta.json" "Crystallography Open Database" "$url" "CC0-1.0" \
    "S. Gražulis et al., Nucleic Acids Res. 40, D420 (2012); при использовании записи цитировать её первоисточник (поле doi в index.tsv)" \
    "Экспериментальные структуры; фильтровать по Robs, celltemp, cellpressure. Записей найдено: $n, CIF: $got (лимит $max_per_material)."
done
