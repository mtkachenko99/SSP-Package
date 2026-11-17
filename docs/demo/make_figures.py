"""Рисунки для демонстрационного документа (данные — результаты ssp_aggregate / ssp_benchmark).

Запуск: software/qc/venv/bin/python docs/demo/make_figures.py
"""
import csv
import glob
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DEMO = os.path.join(ROOT, "software/excitons/examples/demo_c8s3")
RES = os.path.join(DEMO, "results")
OUT = os.path.join(ROOT, "docs/demo/figures")
os.makedirs(OUT, exist_ok=True)

# палитра (dataviz reference, светлая тема): категориальные слоты по порядку; эксперимент — нейтральный
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK, INK2, GRID, EXP = "#1f1f1e", "#5c5b55", "#e6e5e0", "#3d3d3a"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK,
    "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 1.6, "legend.frameon": False,
    "figure.dpi": 200, "savefig.bbox": "tight", "axes.titlesize": 9.5, "axes.titleweight": "bold",
})


def read(path):
    rows = [r for r in csv.reader(l for l in open(path, encoding="utf-8") if not l.startswith("#"))]
    head = rows[0]
    data = []
    for r in rows[1:]:
        try:
            data.append([float(x) if x != "" else np.nan for x in r])
        except ValueError:
            data.append(r)
    return head, data


def cols(path):
    h, d = read(path)
    a = np.array([r for r in d if not isinstance(r[0], str)], dtype=float)
    return {k: a[:, i] for i, k in enumerate(h)}


def xy_file(path):
    x, y = [], []
    for l in open(path, encoding="utf-8"):
        if l.startswith("#") or not l.strip():
            continue
        p = [v for v in l.replace(",", " ").split()]
        try:
            x.append(float(p[0])); y.append(float(p[1]))
        except ValueError:
            pass
    return np.array(x), np.array(y)


def save(fig, name):
    fig.savefig(os.path.join(OUT, name))
    plt.close(fig)
    print("figure", name)


def nm_axis(ax):
    sec = ax.secondary_xaxis("top", functions=(lambda w: 1e7 / np.maximum(w, 1), lambda l: 1e7 / np.maximum(l, 1)))
    sec.set_xlabel("λ, nm", color=INK2)
    sec.tick_params(colors=INK2)


# 1. Калибровка мономера
def fig_calibration():
    d = os.path.join(RES, "1_monomer_calibration")
    a, e = cols(os.path.join(d, "fit_absorption.csv")), cols(os.path.join(d, "fit_emission.csv"))
    fig, (ax, axr) = plt.subplots(2, 1, figsize=(6.3, 4.2), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
    na, ne = a["data"].max(), e["data"].max()
    ax.plot(a["w_cm"], a["data"] / na, "o", ms=2.2, color=EXP, label="absorption (experiment)")
    ax.plot(a["w_cm"], a["model"] / na, color=C[0], label="absorption (model)")
    ax.plot(e["w_cm"], e["data"] / ne, "s", ms=2.0, color="#8a8a84", label="PL (digitized)")
    ax.plot(e["w_cm"], e["model"] / ne, color=C[1], label="PL (model)")
    ax.set_ylabel("normalized intensity")
    ax.legend(loc="upper left", fontsize=8)
    ax.set_title("Step 1. C8S3 monomer in methanol: cumulant-model fit")
    nm_axis(ax)
    axr.plot(a["w_cm"], a["resid"], "o", ms=1.8, color=C[0])
    axr.plot(e["w_cm"], e["resid"], "s", ms=1.8, color=C[1])
    axr.axhline(0, color=INK2, lw=0.8)
    axr.set_ylabel("residual, σ")
    axr.set_xlabel("wavenumber, cm⁻¹")
    axr.set_xlim(15500, 21800)
    save(fig, "fig1_monomer_calibration.png")


# 2. Нанотрубка: уровни модели против эксперимента
TUBES = [("no_dispersion", "no dispersion shift"), ("extended_dipole", "extended dipole + dispersion"),
         ("huckel_charges", "Hückel charges + dispersion"), ("tddft_charges", "TD-DFT charges + dispersion")]


def peak(w, a):
    i = int(np.nanargmax(a))
    return w[i]


def fig_tube():
    lam, ex = xy_file(os.path.join(DEMO, "data/nanotubes_abs_sugar.csv"))
    w_ex = 1e7 / lam
    fig, ax = plt.subplots(figsize=(6.3, 3.6))
    ax.fill_between(w_ex, ex / ex.max(), color="#d9d8d2", lw=0, label="experiment (C8S3 nanotubes)")
    ax.plot(w_ex, ex / ex.max(), color=EXP, lw=1.2)
    rows = []
    for k, (name, lab) in enumerate(TUBES):
        p = os.path.join(RES, "2_nanotube", name, "spectra.csv")
        if not os.path.exists(p):
            continue
        s = cols(p)
        y = s["abs_nm2"] / s["abs_nm2"].max()
        ax.plot(s["w_cm"], y, color=C[k], label=lab)
        rows.append((lab, peak(s["w_cm"], s["abs_nm2"])))
    ax.axvline(19231, color=INK2, lw=0.8, ls=":")
    ax.text(19260, 0.93, "monomer", color=INK2, fontsize=7.5)
    ax.set_xlim(15500, 20500)
    ax.set_ylim(0, 1.08)
    ax.set_xlabel("wavenumber, cm⁻¹")
    ax.set_ylabel("absorption (normalized)")
    ax.set_title("Step 2. Nanotube J-band: prediction without fitting")
    ax.legend(loc="upper left", fontsize=7.5)
    nm_axis(ax)
    save(fig, "fig2_nanotube_levels.png")
    return rows, peak(w_ex, ex)


# 3. Дисперсионные сдвиги сайтов
def fig_site_shifts():
    p = os.path.join(RES, "2_nanotube", "tddft_charges", "site_shifts.csv")
    s = cols(p)
    inner = s["r_xy_nm"] < np.median(s["r_xy_nm"]) * 0 + (s["r_xy_nm"].min() + s["r_xy_nm"].max()) / 2
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.6, 2.9), gridspec_kw={"wspace": 0.35})
    bins = np.linspace(np.percentile(s["total_shift_cm"], 0.5), np.percentile(s["total_shift_cm"], 99.5), 40)
    a1.hist(s["total_shift_cm"][inner], bins=bins, color=C[0], alpha=0.85, label="inner wall")
    a1.hist(s["total_shift_cm"][~inner], bins=bins, color=C[1], alpha=0.75, label="outer wall")
    a1.set_xlabel("site energy shift, cm⁻¹")
    a1.set_ylabel("number of molecules")
    a1.legend(fontsize=7.5)
    a1.set_title("by wall")
    a2.plot(s["z_nm"][inner], s["total_shift_cm"][inner], ".", ms=2, color=C[0])
    a2.plot(s["z_nm"][~inner], s["total_shift_cm"][~inner], ".", ms=2, color=C[1])
    a2.set_xlabel("position along the tube axis, nm")
    a2.set_ylabel("shift, cm⁻¹")
    a2.set_title("along the tube (ends blue-shifted)")
    fig.suptitle("Step 3a. Dispersive site-energy shifts (TD-DFT charges)", fontsize=9.5, fontweight="bold", y=1.06)
    save(fig, "fig3_site_shifts.png")
    return float(np.median(s["total_shift_cm"][inner])), float(np.median(s["total_shift_cm"][~inner]))


# 4. Роль беспорядка
def fig_disorder():
    d = os.path.join(RES, "3_disorder_scan", "tddft_charges")
    a = cols(os.path.join(d, "scan_abs.csv"))
    o = cols(os.path.join(d, "scan_observables.csv"))
    lam, ex = xy_file(os.path.join(DEMO, "data/nanotubes_abs_sugar.csv"))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.6, 2.9), gridspec_kw={"wspace": 0.32})
    keys = [k for k in a if k.startswith("v=")]
    shades = ["#bcd5f4", "#8ab6eb", "#5a96e0", "#2a78d6", "#1d5aa6", "#123d72"]  # одна шкала: σ растёт → темнее
    for k, col in zip(keys, shades):
        a1.plot(a["w_cm"], a[k] / a[k].max(), color=col, lw=1.3, label="σ = " + k[2:])
    a1.plot(1e7 / lam, ex / ex.max(), color=EXP, lw=1.0, ls="--", label="experiment")
    a1.set_xlim(16000, 18500)
    a1.set_xlabel("wavenumber, cm⁻¹"); a1.set_ylabel("absorption (normalized)")
    a1.set_title("spectra for σ (cm⁻¹)")
    a1.legend(fontsize=6.8, ncol=1)
    a2.plot(o["value"], o["abs_fwhm_cm"], "o-", color=C[0], ms=4, label="model")
    a2.axhline(423, color=EXP, ls="--", lw=1)
    a2.text(5, 440, "measured: 423 cm⁻¹", color=INK2, fontsize=7.5)
    a2.axvline(124.5, color=C[1], ls=":", lw=1)
    a2.text(130, 880, "monomer\nσ_inh", color=C[1], fontsize=7.5)
    a2.set_xlabel("static disorder σ, cm⁻¹"); a2.set_ylabel("J-band FWHM, cm⁻¹")
    a2.set_title("band width")
    fig.suptitle("Step 3b. Disorder controls the J-band width", fontsize=9.5, fontweight="bold", y=1.04)
    save(fig, "fig4_disorder.png")
    return list(zip(o["value"], o["abs_fwhm_cm"], o["abs_peak_cm"]))


# 6. Распад у золотой сферы: модель против точной теории Ми
def fig_decay():
    h, d = read(os.path.join(ROOT, "data/results/benchmarks/purcell_vs_mie.csv"))
    conv = lambda v: v if isinstance(v, float) else (float(v) if v.replace(".", "", 1).replace("-", "", 1).replace("e", "", 1).isdigit() else v)
    rows = [{k: conv(v) for k, v in zip(h, r)} for r in d]
    fig, axs = plt.subplots(1, 2, figsize=(6.6, 2.9), gridspec_kw={"wspace": 0.32})
    for ax, orient, ttl in [(axs[0], "radial", "dipole ⊥ surface"), (axs[1], "tangential", "dipole ∥ surface")]:
        for k, R in enumerate([5.0, 10.0, 20.0]):
            sel = [r for r in rows if r["radius_nm"] == R and r["orientation"] == orient and r["lambda_nm"] == 650.0]
            g = np.array([r["gap_nm"] for r in sel])
            ax.plot(g, [r["exact_qy"] for r in sel], "o", color=C[k], ms=4, mfc="none")
            ax.plot(g, [r["model_qy"] for r in sel], "-", color=C[k], lw=1.4, label=f"R = {R:g} nm")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("gap to the surface, nm"); ax.set_title(ttl)
    axs[0].set_ylabel("quantum yield (intrinsic = 1)")
    axs[0].legend(fontsize=7.5, title="line: model, ○: exact Mie", title_fontsize=7)
    fig.suptitle("Molecule near a gold sphere in water, λ = 650 nm: quenching and emission", fontsize=9.5, fontweight="bold", y=1.04)
    save(fig, "fig6_decay_vs_mie.png")


# 5. Гибрид: слой C8S3 на Au-сфере
def fig_hybrid():
    d = os.path.join(RES, "4_hybrid_scan")
    m = cols(os.path.join(d, "shell", "scan_abs_diff.csv"))
    keys = [k for k in m if k.startswith("v=")]
    vals = np.array([float(k[2:]) for k in keys])
    Z = np.array([m[k] for k in keys])
    br = cols(os.path.join(d, "shell", "scan_branches.csv"))
    qy, qy0 = [], []
    for i in range(len(vals)):
        line = open(os.path.join(d, f"scan_{i:02d}", "shell", "hybrid_emission.csv"), encoding="utf-8").read().splitlines()[1]
        a = line.split("QY")[1].split(";")[0].split("vs")
        qy.append(float(a[0])); qy0.append(float(a[1]))
    h, cp = read(os.path.join(d, "shell", "scan_coupling.csv"))
    coup = {k: v for k, v in zip(h, cp[0])}
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.6, 3.0), gridspec_kw={"wspace": 0.62, "width_ratios": [1.25, 1]})
    lim = np.nanmax(np.abs(Z))
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("div", ["#1d5aa6", "#9cc0ea", "#f2f1ec", "#f4a77f", "#b9461a"])
    im = a1.pcolormesh(m["w_cm"], vals, Z, cmap=cmap, vmin=-lim, vmax=lim, shading="auto")
    a1.plot(br["w_p"], vals, color=INK, lw=1, ls="--", label="plasmon")
    a1.plot(vals, vals, color=INK2, lw=0.8, ls=":", label="exciton")
    a1.set_xlim(15500, 22500)
    a1.set_xlabel("wavenumber, cm⁻¹"); a1.set_ylabel("exciton energy (scan), cm⁻¹")
    a1.set_title("aggregate contribution to absorption")
    cb = fig.colorbar(im, ax=a1, pad=0.02); cb.set_label("Δσ, nm²", color=INK2)
    a1.legend(fontsize=7, loc="lower right")
    a2.semilogy(vals, qy0, "o-", color=C[0], ms=3.5, label="aggregate alone")
    a2.semilogy(vals, qy, "s-", color=C[1], ms=3.5, label="1 nm from gold")
    a2.set_xlabel("exciton energy, cm⁻¹"); a2.set_ylabel("quantum yield")
    a2.set_title("emission quenching")
    a2.legend(fontsize=7)
    fig.suptitle("Step 4. C8S3-like layer on a gold sphere (R = 20 nm, water)", fontsize=9.5, fontweight="bold", y=1.04)
    save(fig, "fig5_hybrid_map.png")
    return dict(g=float(coup["g_cm"]), g_hi=float(coup["g_hi_cm"]), gamma_p=float(coup["gamma_p_cm"]), gamma_x=float(coup["gamma_x_cm"]),
                C=float(coup["cooperativity"]), omega_p=float(np.nanmean(br["w_p"])), qy_min=min(qy), qy_max=max(qy), qy0=float(np.mean(qy0)))


# 7. HEOM против Редфилда
def fig_heom():
    p = os.path.join(ROOT, "software/excitons/examples/results/heom_benchmark/chain6/heom_populations.csv")
    s = cols(p)
    fig, ax = plt.subplots(figsize=(5.6, 2.9))
    show = [(5, "highest state (initial)"), (0, "lowest state")]
    for k, (i, lab) in enumerate(show):
        ax.plot(s["t_ps"], s[f"heom_{i}"], color=C[k], lw=1.8, label=lab + ": HEOM")
        ax.plot(s["t_ps"], s[f"redfield_{i}"], color=C[k], lw=1.2, ls="--", label=lab + ": Redfield")
    ax.set_xlabel("time, ps"); ax.set_ylabel("population")
    ax.set_title("Exciton relaxation in a 6-molecule chain: exact HEOM vs Redfield")
    ax.legend(fontsize=7, ncol=2)
    save(fig, "fig7_heom.png")


# 8. Формы агрегатов
def fig_shapes():
    base = os.path.join(ROOT, "software/excitons/examples/results/notebook_demo")
    names = sorted(n for n in os.listdir(base) if os.path.isdir(os.path.join(base, n)))
    fig, ax = plt.subplots(figsize=(5.6, 2.9))
    for k, n in enumerate([n for n in names if n != "double_wall_tube"][:4]):  # полоса трубки выходит за окно примера
        sp = cols(os.path.join(base, n, "spectra.csv"))
        ax.plot(sp["w_cm"], sp["abs_nm2"] / sp["abs_nm2"].max(), color=C[k], label=n.replace("_", " "))
    ax.set_xlabel("wavenumber, cm⁻¹"); ax.set_ylabel("absorption (normalized)")
    ax.set_title("Absorption of aggregates of different shape")
    ax.legend(fontsize=7.5)
    nm_axis(ax)
    save(fig, "fig8_shapes.png")


def calib_numbers():
    h, d = read(os.path.join(RES, "1_monomer_calibration", "fit_parameters.csv"))
    out = {}
    for r in d:
        if isinstance(r[0], str):
            name = r[0]
            try:
                out[name] = float(r[1]); out[name + "_err"] = float(r[2]) if r[2] not in ("", None) else 0.0
            except (ValueError, IndexError):
                pass
    chi = [l for l in open(os.path.join(RES, "1_monomer_calibration", "fit_parameters.csv")) if "chi2_red" in l][0]
    out["chi2"] = float(chi.split("=")[-1])
    out["stokes"] = out.get("stokes_peak_cm (model)", 0.0)
    return out


if __name__ == "__main__":
    import json
    fig_calibration()
    tube, tube_exp = fig_tube()
    sin, sout = fig_site_shifts()
    dis = fig_disorder()
    fig_decay()
    hyb = fig_hybrid()
    fig_heom()
    fig_shapes()
    hyb["text"] = (f"Сила связи слоя с дипольным плазмоном g = {hyb['g']:.0f} см⁻¹ (верхняя граница интервала {hyb['g_hi']:.0f} см⁻¹) при ширинах "
                   f"плазмона {hyb['gamma_p']:.0f} и экситонной полосы {hyb['gamma_x']:.0f} см⁻¹: кооперативность C = {hyb['C']:.3f}, "
                   f"режим слабой связи — гибридные ветви не разрешаются, на карте видны лишь интерференционные провал и подъём "
                   f"у резонанса. Квантовый выход падает с {hyb['qy0']:.2f} до {hyb['qy_min']:.1e}–{hyb['qy_max']:.1e}: на расстоянии 1 нм "
                   f"преобладает тушение металлом. Для сильной связи нужен плотный упорядоченный J-агрегат (узкая полоса и "
                   f"большая плотность диполей) — это можно исследовать той же программой (пример dense_ligand_shell).")
    json.dump({"cal": calib_numbers(), "tube": [{"label": l, "peak": float(p)} for l, p in tube], "tube_exp": float(tube_exp),
               "shift_in": sin, "shift_out": sout, "disorder": [[float(a), float(b), float(c)] for a, b, c in dis],
               "hybrid": hyb}, open(os.path.join(os.path.dirname(__file__), "numbers.json"), "w"), ensure_ascii=False, indent=1)
    print("numbers.json written")
