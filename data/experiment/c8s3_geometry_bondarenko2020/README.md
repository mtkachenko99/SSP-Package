# C8S3: атомная геометрия (DFT) и заряды

- **Источник:** A. S. Bondarenko, I. Patmanidis, R. Alessandri, P. C. T. Souza, A. H. de Vries, S. J. Marrink,
  J. Knoester, «C8S3 molecular geometry and charge distributions», figshare (2020),
  doi:10.6084/m9.figshare.12369335, лицензия **CC BY 4.0**; к работе Bondarenko et al., ChemRxiv 2020,
  doi:10.26434/chemrxiv.12279956.
- **Скачано:** 2026-11-15, `raw/`: `gs_B3LYP-631Gdp_optimized.xyz` (108 атомов, sha256 15e8494b110a2337…), CHELPG-заряды основного и
  возбуждённого состояний (ωB97XD/6-31G(d,p)) — это **статические**, не переходные заряды.
- **processed/c8s3_huckel_tresp.txt:** переходные заряды Хюккеля (`ssp_aggregate --huckel-tresp … 1`): 21 π-атом,
  24 π-электрона, щель HOMO–LUMO 0.82|β|.
- **tddft_camb3lyp_631gs/:** хромофор с метилами вместо цепей (48 атомов), PySCF CAM-B3LYP/6-31G*, TDA, переходные
  заряды по ЭСП переходной плотности (`software/qc/c8s3_tresp.py`).
