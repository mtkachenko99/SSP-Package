#!/usr/bin/env bash
# Выгрузка снапшота JARVIS-DFT 3D (NIST, Figshare) и извлечение записей для materials.tsv.
# Полный архив (~41 МБ zip) хранится в кэше; распакованный JSON не сохраняется.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
need curl jq unzip

snapshot="jdft_3d-12-12-2022"
url="https://ndownloader.figshare.com/files/38521619"
dir="$SSP_REF/jarvis_dft_3d"
zip="$dir/cache/$snapshot.json.zip"
mkdir -p "$dir/cache"

if [[ ! -s "$zip" ]]; then
  log "JARVIS: скачивание $snapshot"
  curl -sSfL --retry 3 -A "$SSP_UA" "$url" -o "$zip.part" && mv "$zip.part" "$zip"
fi

formulas="$(materials | cut -f3 | jq -Rn '[inputs]')"
log "JARVIS: извлечение записей для $(jq length <<<"$formulas") формул"
unzip -p "$zip" | jq --argjson F "$formulas" '[ .[] | select(.formula as $f | any($F[]; . == $f)) ]' \
  > "$dir/starter_set.json"

# Краткая сводка: jid, формула, группа, щели, эфф. массы, модули, энергия над выпуклой оболочкой
jq -r '["jid","formula","spg_number","spg_symbol","optb88vdw_bandgap","mbj_bandgap","hse_gap","avg_elec_mass","avg_hole_mass","bulk_modulus_kv","shear_modulus_gv","epsx","mepsx","ehull","exfoliation_energy"],
       (.[] | [.jid,.formula,.spg_number,.spg_symbol,.optb88vdw_bandgap,.mbj_bandgap,.hse_gap,.avg_elec_mass,.avg_hole_mass,.bulk_modulus_kv,.shear_modulus_gv,.epsx,.mepsx,.ehull,.exfoliation_energy])
       | map(if . == null then "" else tostring end) | @tsv' "$dir/starter_set.json" > "$dir/starter_set.tsv"
log "JARVIS: записей $(jq length "$dir/starter_set.json") → $dir/starter_set.tsv"

write_meta "$dir.meta.json" "JARVIS-DFT 3D (NIST)" "$url" \
  "NIST public data (проверить условия на jarvis.nist.gov)" \
  "K. Choudhary et al., npj Comput. Mater. 6, 173 (2020)" \
  "Полные записи (включая частотные ε(ω), упругие тензоры, ИК/КР, где есть) — starter_set.json; сводка — starter_set.tsv. 'na' = не рассчитано." \
  "OptB88vdW (структуры, ε(ω) в RPA без локальных полей), TBmBJ (mbj_*), без экситонных эффектов" "$snapshot"
