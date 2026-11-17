# software/ — код проекта

C++20 (расчётные ядра), Qt6 (GUI, сетевой доступ к базам), bash (сборка, выгрузка данных, пакетные расчёты).
Сборка — CMake-суперпроект: модуль подключается автоматически, когда в его каталоге появляется `CMakeLists.txt`.

```bash
scripts/build.sh            # Release
scripts/build.sh Debug -DSSP_BUILD_GUI=OFF
```

## Модули

| Каталог | Назначение | Методология / обзор |
|---|---|---|
| `common/` | Библиотека `ssp_common`: единицы и константы, сетки, ε(ω)-материалы (таблицы RII, Друде–Лоренц, Крейбиг, LST, лоренцевы экситоны), структуры (CIF через gemmi, симметрия через spglib), ввод/вывод HDF5/CSV/JSON, provenance-метаданные | M2, M3.1, M9 |
| `data_fetch/` | Выгрузка справочных данных: OPTIMADE, COD, refractiveindex.info, JARVIS, MP, Raman DB, RRUFF, C2DB (bash + Qt Network) | M2.1–M2.2, [01_databases](../docs/01_databases.md) |
| `plasmonics/` | Квазистатика/Ганс/MLWA, многослойный Ми (+HDM, d-параметры, QNM-полюса, диполь у сферы), позже T-матрица, DDA; обёртки Meep/ADDA/SCUFF-EM | [plasmonics_polaritons](../docs/methods/plasmonics_polaritons.md) |
| `excitons/` | EMA (Брус, непараболическая, диэлектрическое ограничение, core/shell), Эллиотт, Келдыш–Рытова, тонкая структура, связанные осцилляторы, Тэвис–Каммингс/Линдблад | [excitons_phonons_spectra](../docs/methods/excitons_phonons_spectra.md) §1, plasmonics §5 |
| `spectroscopy/absorption/` | Ансамблевые спектры поглощения/экстинкции (P(R), ориентации, аппаратная функция) | M3.2 |
| `spectroscopy/fluorescence/` | Хуанг–Рис (производящая функция), Сегалл Γ(T), кинетика, усиление/тушение у НЧ | excitons §2–3, plasmonics §6 |
| `spectroscopy/raman/` | Правила отбора, PCM, BPM, A-член Альбрехта, SERS \|E\|⁴, чтение Phonopy | excitons §4 |
| `spectroscopy/brillouin/` | Лэмб (свободная/в матрице/core–shell), Дюваль, Рэлей–Ритц (Visscher), фотоупругие сечения | excitons §5 |
| `validation/` | Подгонка (LM, CMA-ES), MCMC, модель расхождения (GP), анализ остатков (PCA/NMF, автокорреляция), Соболь, обращение P(R) | M4, M6–M8 |
| `gui/` | Qt6-приложение: проекты, запуск расчётов, графики «модель + остатки», карты | [04_results_presentation](../docs/04_results_presentation.md) |
| `scripts/` | bash: сборка, пакетные прогоны, выгрузка данных | — |
| `third_party/` | Внешние библиотеки (submodule/vendored), если нет в системе | — |
| `cmake/` | CMake-модули поиска зависимостей | — |

## Зависимости

Установлено в системе (проверено 2026-11-11): GCC 16, CMake, Ninja, Qt 6.12 (Core, Widgets), FFTW3, HDF5.

Нужно доустановить (Arch: `pacman`/AUR) по мере появления модулей:

| Библиотека | Для чего | Пакет (Arch) |
|---|---|---|
| Eigen 3.4/5.x | Линейная алгебра, задачи на СЗ (установлен) | `eigen` |
| Boost.Math | Сферические функции Бесселя/Ханкеля, квадратуры, поиск корней | `boost` |
| gemmi | Парсер CIF | AUR / vendored (header-only) |
| spglib | Симметрия кристаллов | `spglib` |
| yaml-cpp | YAML refractiveindex.info | `yaml-cpp` |
| nlohmann-json | JSON (вне Qt-кода) | `nlohmann-json` |
| Ceres Solver | Нелинейный МНК | `ceres-solver` |
| GoogleTest / Catch2 | Тесты | `gtest` / `catch2` |
| Qt Graphs или QCustomPlot | Графики в GUI | `qt6-graphs` / vendored |

## Соглашения

- Пространство имён `ssp::<module>`; заголовки `common/include/ssp/...`.
- Внутренние единицы: СИ, энергия — эВ на интерфейсах; явные суффиксы в именах (`lambda_nm`, `E_eV`).
- Каждый исполняемый расчёт принимает конфиг (JSON) и пишет HDF5 с `/config` и provenance (M9).
- Каждый физический модуль имеет тесты аналитических пределов и бенчмарков (M4).

## Выгрузка справочных данных

```bash
software/data_fetch/fetch_all.sh          # refractiveindex.info + JARVIS + COD для materials.tsv
COD_MAX=50 software/data_fetch/fetch_cod.sh  # ограничить число CIF на материал
```
Список материалов — `data_fetch/materials.tsv` (TAB-разделённый; формулы COD в нотации Хилла через пробел).
Требуются `curl`, `jq`, `git`, `unzip`. Скрипты идемпотентны: уже скачанное не перекачивается.

## ssp_aggregate — агрегаты произвольной формы (+ плазмонная НЧ)

Сборка и тесты:
```bash
software/scripts/build.sh
ctest --test-dir software/build
```
Запуск (пути в конфиге, включая output_dir, — относительно файла конфига):
```bash
cd software/excitons/examples && ../../build/excitons/ssp_aggregate shapes_with_gold.json
software/build/excitons/ssp_aggregate --list-geometries
```
Примеры: `notebook_demo.json` (формы и параметры testcalc.ipynb, связь legacy_J0), `shapes_with_gold.json`
(4 формы + Au-сфера 15 нм), `heom_benchmark.json` (кумулянт+Редфилд vs HEOM), `shell_scan.json` (скан по оболочке),
`dense_ligand_shell.json` (плотная доменная оболочка на Au с лигандами, этап C). Модель и проверки —
[docs/methods/aggregate_plexciton_model.md](../docs/methods/aggregate_plexciton_model.md).

Выход (на каждую форму): `geometry.csv`, `spectra.csv` (σ_abs и ФЛ ± SE), `observables.csv` (17 параметров),
`states_representative.csv`, `dynamics.csv`, при НЧ — `plexciton.csv`, `molecules_near_field.csv`,
`hybrid_emission.csv`, `np_mode.csv`, при модели Хольстейна — `vibronic_absorption.csv`, при HEOM —
`heom_absorption.csv`, `heom_populations.csv`; общая `summary.csv` для сравнения форм.
Eigen: поддерживаются 3.4 и 5.x.

## ssp_gui — графическая оболочка (Qt6 Widgets, интерфейс на английском)

```bash
software/build/gui/ssp_gui                          # открывает пример shapes_with_gold.json
software/build/gui/ssp_gui path/to/config.json
```
- **Слева — параметры**: «Формы» (список агрегатов, тип, параметры, размещение, 3D-превью с наночастицей; вращение
  мышью), «Модель», «Расчёт», «Наночастица», «HEOM», «JSON» (полный конфиг; форма ↔ JSON синхронны). Конфиг —
  тот же JSON, что читает `ssp_aggregate`; «Сохранить» пишет его на диск.
- **▶ Запустить (F5)** — расчёт в фоне, журнал и прогресс; «Остановить» — после текущей формы.
- **Справа — результаты**: спектры всех форм (±SE, нормировка, шкала λ), «Агрегат + НЧ» (σ_ext/abs/sca системы и
  НЧ, вклад агрегата), кинетика, HEOM-эталон с панелью невязки, ближнее поле (усиление, тушение, 3D-раскраска),
  сводная таблица форм, журнал. Графики: выделение мышью — увеличение, колесо — масштаб, двойной щелчок — весь
  диапазон, правая кнопка — экспорт PNG/SVG/CSV.
- **Сравнение с экспериментом** (вкладка «Спектры»): загрузка таблицы x, y[, σ] (см⁻¹, нм или эВ), подгонка
  масштаба, базовой линии и сдвига δ; панель остатков (в σ), χ²_ν, Дарбин–Уотсон — признак неучтённых эффектов (M7).
- Для автоматизации/проверки: `--results DIR`, `--experiment FILE --exp-unit nm`, `--fit offset,shift`, `--run`,
  `--tab N`, `--screenshot out.png` (с `QT_QPA_PLATFORM=offscreen`).

### Оболочка на наночастице и сканы
Форма `particle_shell` — слой молекул на поверхности НЧ (полуоси берутся из раздела `nanoparticle`). Раздел конфига
```json
"scan": {"enabled": true, "parameter": "exciton.monomer_energy_cm", "from": 16500, "to": 21500, "steps": 11,
         "also_set": "nanoparticle.semi_axis_z_nm"}
```
запускает расчёт для каждого значения (`output_dir/scan_NN`) и пишет в `output_dir/<форма>/`: карты `scan_*.csv`,
`scan_branches.csv` (пики плазмона, экситона, гибридные пики, ветви модели), `scan_coupling.csv` (g, 2g, критерии
сильной связи для поглощения/экстинкции/рассеяния), `scan_observables.csv`. Пример: `examples/shell_scan.json`.
В GUI: вкладка параметров **Scan**, вкладка результатов **Scan** (тепловая карта с ветвями, таблица критериев,
наблюдаемая от параметра).

### Этап C: реалистичные оболочки, лиганды, излучение (v0.5)
Новые ключи конфига (все необязательные; полный пример — `examples/dense_ligand_shell.json`, N ≈ 1000, ~5–8 мин на 16 потоках):
```json
"configurations": [{"geometry_type": "particle_shell",
   "geometry_params": {"gap": 0.5, "density_per_nm2": 0.6, "dipole_mode": 1, "domains": 12, "domain_jitter_deg": 10}}],
"exciton": {"disorder_corr_length_nm": 2.0, "orientation_jitter_deg": 5, "position_jitter_nm": 0.05,
            "relaxation": "modified_redfield",
            "coupling": {"model": "transition_charges", "charges_file": "dye_tresp.txt"}},
"holstein": {"enabled": true, "omega_cm": 1300, "huang_rhys": 0.6, "nu_max": 4, "realizations": 20},
"nanoparticle": {"shell_thickness_nm": 1.0, "shell_n": 1.45, "d_perp_nm": 0.0, "d_perp_im_nm": 0.0,
                 "fast_l": "auto", "size_sigma_rel": 0.1, "size_nodes": 5, "emission": true}
```
- `particle_shell` ставит молекулы над **внешней** поверхностью (ядро + `shell_thickness_nm`); `density_per_nm2`
  вместо `N`; `domains` — упорядоченные домены (диполи по направлению домена с разбросом `domain_jitter_deg`).
- `charges_file` — столбцы `x_nm y_nm z_nm q` в системе молекулы (x вдоль диполя), нормируются на μ.
- `relaxation: modified_redfield` — для N ≲ 300 (стоимость O(N²·n_t)).
- `fast_l`: `"auto"` (по умолчанию: быстрый решатель при N ≥ 600), `0` — точный N×N, `L` — мультиполи l ≤ L точно.
- Новые файлы: `hybrid_emission.csv` (ФЛ у НЧ и без неё, QY, время жизни), `vibronic_absorption.csv`;
  в `molecules_near_field.csv` — `rad_enhancement`, `nonrad_to_np_ps`, `quantum_yield`.
- GUI: вкладка «Aggregate + NP» — выбор «PL: hybrid vs aggregate without NP»; «Near field» — тушение / Γ_rad/Γ₀ /
  квантовый выход по зазору; на «Spectra» пунктиром — поглощение в модели Хольстейна.

### Этап B: валидация и калибровка (v0.6)
```bash
software/build/validation/ssp_benchmark purcell <Au.yml|drude> out.csv      # модель vs точная теория Ми
software/build/validation/ssp_benchmark extinction <Au.yml>... out.csv     # λmax и ε Au-сфер vs паспортные данные
software/build/excitons/ssp_aggregate excitons/examples/monomer_calibration.json   # калибровка по мономеру
```
Калибровочный JSON — раздел `monomer_fit`: `absorption_file` (x, y[, σ]; `absorption_unit` nm|cm-1|eV,
`absolute_extinction` — y в М⁻¹см⁻¹ → μ и τ_rad), `emission_file` (`emission_per_wavelength` — пересчёт с λ²),
`n_medium`, `temperature_K`, `wmin_cm`/`wmax_cm`, `start` (E0_cm, lambda_cm, cutoff_cm, sigma_inh_cm, modes),
`fit` (флаги E0, lambda, cutoff, sigma_inh, mode_S, mode_omega, moment_start), `output_dir`. Выход:
`fit_parameters.csv`, `fit_absorption.csv`, `fit_emission.csv`, `aggregate_config_snippet.json`.
В GUI: **File → Calibrate monomer…**, вкладка **Calibration** (данные, модель, остатки, параметры, кнопка
«Apply to the aggregate model»). Данные в `examples/monomer_synthetic/` — синтетические (сгенерированы моделью).

### Дисперсионный сдвиг энергий сайтов (v0.7)
```json
"exciton": {"dispersion": {"enabled": true, "F_eV": 1.43, "n_aggregate": 1.78, "n_solvent": 1.329,
            "dipole_length_nm": 1.0, "softening_nm": 0.3, "cutoff_nm": 3.0, "extra_shift_cm": 0}}
```
Энергия сайта = энергия мономера в растворителе + возврат к газовой фазе + Q·Σ_k V_mk² (Megow et al. 2015).
Геометрия `brickwork_tube` — трубка из кирпичной решётки с бесшовной свёрткой (параметры C8S3 по умолчанию).
Пример с данными: `data/experiment/doria2022_c8s3_nanotubes/aggregate_prediction.json` (и `_nodisp`).
