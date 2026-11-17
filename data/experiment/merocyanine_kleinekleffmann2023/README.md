# Мероцианин MC1: мономер в 5 растворителях, H- и J-агрегаты

- **Источник:** L. Kleine-Kleffmann, A. Schulz, V. Stepanenko, F. Würthner, «Growth of Merocyanine Dye J-Aggregate
  Nanosheets by Living Supramolecular Polymerization», *Angew. Chem.* (2023), doi:10.1002/anie.202314667.
- **Данные:** Zenodo, doi:10.5281/zenodo.10066797, лицензия **CC BY 4.0**.
- **Скачано:** 2026-11-15; `raw/AdditionalData_anie.202314667.xlsx` (0.89 МБ). Чтение — stdlib-скриптом (без openpyxl).
- **processed/**: лист «Figure 2 + Figure S1» (293 K, c = 10⁻⁵ M) — мономер в 1,4-диоксане, хлороформе, THF, MeCN,
  DMSO и агрегаты Agg1 (H, 499 нм) и Agg2 (J, 625 нм) в MeCN/вода 1:9. Столбец подписан «molar extinction
  coefficient (a.u.)», но значения (5.8–6.6·10⁴ в максимуме) — молярный коэффициент экстинкции, М⁻¹см⁻¹.
- **Ограничения:** спектров ФЛ нет; краситель сильно сольватохромен (максимум 513 нм в диоксане — 569 нм в DMSO),
  поэтому энергия сайта в водном агрегате известна с неопределённостью ~±1000 см⁻¹.
- **Калибровка:** `calibration_monomer_mecn.json`, `calibration_monomer_dioxane.json` →
  `data/results/calibration_merocyanine_*`.
