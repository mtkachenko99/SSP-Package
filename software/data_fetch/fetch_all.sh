#!/usr/bin/env bash
# Выгрузка всех открытых справочных данных для стартового набора материалов.
# Источники с ручным доступом (RRUFF/AMCSD — CAPTCHA) сюда не входят, см. docs/01_databases.md.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
"$here/fetch_refractiveindex.sh"
"$here/fetch_jarvis.sh"
"$here/fetch_cod.sh"
