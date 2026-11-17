"""TD-DFT переходные заряды (TrEsp) хромофора C8S3 по геометрии Bondarenko et al. 2020.

Цепи у атомов N (октил, сульфопропил) заменяются метилами: первый sp3-углерод у N сохраняется, остальное
отбрасывается, разорванные связи этого углерода закрываются H (1.09 Å вдоль прежней связи). Затем
CAM-B3LYP/6-31G* (газовая фаза), TDA, несколько синглетов; переходная плотность S0→Sn подгоняется точечными
зарядами на тяжёлых атомах π-системы и Cl по электростатическому потенциалу (оболочки 1.4–2.0 ван-дер-ваальсовых
радиусов, ограничение Σq = 0, метод Лагранжа). Вывод: x y z q (нм) в системе исходной геометрии и сводка.

Использование: venv/bin/python c8s3_tresp.py geometry.xyz out_dir [xc] [basis] [nstates]
"""
import json
import os
import sys
import time

import numpy as np
from pyscf import dft, gto, lib, tdscf

ANG = 1.0
VDW = {"H": 1.2, "C": 1.7, "N": 1.55, "O": 1.52, "S": 1.8, "Cl": 1.75}


def read_xyz(path):
    lines = open(path).read().splitlines()
    atoms = []
    for line in lines[2:]:
        p = line.split()
        if len(p) >= 4:
            atoms.append((p[0], np.array(list(map(float, p[1:4])))))
    return atoms


def neighbours(atoms):
    n = len(atoms)
    nb = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            hv = atoms[i][0] != "H" and atoms[j][0] != "H"
            if np.linalg.norm(atoms[i][1] - atoms[j][1]) < (1.9 if hv else 1.25):
                nb[i].append(j)
                nb[j].append(i)
    return nb


def cap_chromophore(atoms):
    nb = neighbours(atoms)
    keep = set()
    sp2 = {i for i, a in enumerate(atoms) if a[0] == "C" and len(nb[i]) == 3}
    pi = set(sp2) | {i for i, a in enumerate(atoms) if a[0] == "N" and any(j in sp2 for j in nb[i])}
    # π-система + Cl на кольцах + H на π-углеродах
    for i in pi:
        keep.add(i)
        for j in nb[i]:
            if atoms[j][0] in ("H", "Cl"):
                keep.add(j)
    new = [atoms[i] for i in sorted(keep)]
    # метилы: первый sp3-углерод у каждого N, его H сохраняются, тяжёлые соседи → H
    for i in pi:
        if atoms[i][0] != "N":
            continue
        for c in nb[i]:
            if c in pi or atoms[c][0] != "C":
                continue
            new.append(atoms[c])
            for k in nb[c]:
                if k == i:
                    continue
                if atoms[k][0] == "H":
                    new.append(atoms[k])
                else:
                    v = atoms[k][1] - atoms[c][1]
                    new.append(("H", atoms[c][1] + 1.09 * v / np.linalg.norm(v)))
    return new, len(pi)


def esp_grid(mol, layers=(1.4, 1.6, 1.8, 2.0), density=1.0):
    """Точки на сферах вокруг атомов (Merz–Kollman), без точек внутри других сфер."""
    coords = mol.atom_coords(unit="Angstrom")
    elems = [mol.atom_symbol(i) for i in range(mol.natm)]
    pts = []
    for f in layers:
        for i, (r0, e) in enumerate(zip(coords, elems)):
            R = f * VDW.get(e, 1.7)
            n = max(10, int(4 * np.pi * R * R * density))
            k = np.arange(n) + 0.5
            phi = np.arccos(1 - 2 * k / n)
            th = np.pi * (1 + 5 ** 0.5) * k
            s = r0 + R * np.stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)], 1)
            ok = np.ones(n, bool)
            for j, (r1, e1) in enumerate(zip(coords, elems)):
                if j != i:
                    ok &= np.linalg.norm(s - r1, axis=1) > f * VDW.get(e1, 1.7) - 1e-6
            pts.append(s[ok])
    return np.concatenate(pts)


def transition_potential(mol, T, grid_ang, batch=200):
    """V(r) = −Σ T_μν ⟨μ|1/|r − r'||ν⟩ (атомные единицы), T — переходная матрица плотности в АО."""
    g = grid_ang / 0.52917721092
    V = np.empty(len(g))
    for s in range(0, len(g), batch):
        ints = mol.intor("int1e_grids", grids=g[s:s + batch])
        V[s:s + batch] = -np.einsum("gij,ij->g", ints, T)
    return V, g


def fit_charges(centres_bohr, grid_bohr, V):
    A = 1.0 / np.linalg.norm(grid_bohr[:, None, :] - centres_bohr[None, :, :], axis=2)
    n = centres_bohr.shape[0]
    M = np.zeros((n + 1, n + 1))
    M[:n, :n] = A.T @ A
    M[:n, n] = M[n, :n] = 1.0
    rhs = np.concatenate([A.T @ V, [0.0]])
    q = np.linalg.solve(M, rhs)[:n]
    rrms = np.sqrt(np.mean((A @ q - V) ** 2)) / np.sqrt(np.mean(V ** 2))
    return q, rrms


def main():
    xyz, out = sys.argv[1], sys.argv[2]
    xc = sys.argv[3] if len(sys.argv) > 3 else "camb3lyp"
    basis = sys.argv[4] if len(sys.argv) > 4 else "6-31g*"
    nstates = int(sys.argv[5]) if len(sys.argv) > 5 else 4
    os.makedirs(out, exist_ok=True)
    lib.num_threads(int(os.environ.get("OMP_NUM_THREADS", "16")))
    atoms, npi = cap_chromophore(read_xyz(xyz))
    with open(os.path.join(out, "chromophore_capped.xyz"), "w") as f:
        f.write(f"{len(atoms)}\nC8S3 chromophore, alkyl/sulfopropyl chains replaced by methyl\n")
        for e, r in atoms:
            f.write(f"{e} {r[0]:.6f} {r[1]:.6f} {r[2]:.6f}\n")
    mol = gto.M(atom=[(e, r) for e, r in atoms], basis=basis, charge=1, spin=0, verbose=4,
                output=os.path.join(out, "pyscf.log"))
    t0 = time.time()
    mf = dft.RKS(mol, xc=xc)
    mf.conv_tol = 1e-8
    mf.kernel()
    td = tdscf.TDA(mf)
    td.nstates = nstates
    td.kernel()
    t1 = time.time()
    mo = mf.mo_coeff
    occ = mf.mo_occ > 0
    heavy = [i for i in range(mol.natm) if mol.atom_symbol(i) != "H"]
    grid = esp_grid(mol)
    res = {"xc": xc, "basis": basis, "natoms": mol.natm, "n_pi": npi, "scf_energy_Eh": float(mf.e_tot),
           "minutes": (t1 - t0) / 60, "states": []}
    dips = td.transition_dipole()
    for s in range(nstates):
        X = td.xy[s][0]  # nocc × nvir
        T = np.sqrt(2.0) * mo[:, occ] @ X @ mo[:, ~occ].T  # синглет, спин-адаптированная
        T = 0.5 * (T + T.T)
        V, gb = transition_potential(mol, T, grid)
        cb = mol.atom_coords()[heavy]
        q, rrms = fit_charges(cb, gb, V)
        mu_fit = (q[:, None] * cb).sum(0) * 2.541746  # e·бор → Д
        mu_td = np.array(dips[s]) * 2.541746
        e_ev = float(td.e[s] * 27.211386)
        res["states"].append({"state": s + 1, "E_eV": e_ev, "E_cm": e_ev * 8065.544, "mu_TD_D": float(np.linalg.norm(mu_td)),
                              "mu_fit_D": float(np.linalg.norm(mu_fit)), "esp_rel_rms": float(rrms)})
        with open(os.path.join(out, f"tresp_S{s + 1}.txt"), "w") as f:
            f.write(f"# TrEsp S0->S{s + 1}: {xc}/{basis}, TDA, gas phase; E = {e_ev:.4f} eV; |mu_TD| = "
                    f"{np.linalg.norm(mu_td):.3f} D, |mu_fit| = {np.linalg.norm(mu_fit):.3f} D; ESP rel. rms {rrms:.3f}\n")
            f.write("# geometry: Bondarenko et al. 2020 (figshare 10.6084/m9.figshare.12369335, CC BY 4.0), chains -> methyl\n")
            f.write("# x_nm y_nm z_nm q_e (frame of the source geometry)\n")
            for k, i in enumerate(heavy):
                r = mol.atom_coord(i, unit="Angstrom") * 0.1
                f.write(f"{r[0]:.6f} {r[1]:.6f} {r[2]:.6f} {q[k]:.6f}\n")
    json.dump(res, open(os.path.join(out, "summary.json"), "w"), indent=2)
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
