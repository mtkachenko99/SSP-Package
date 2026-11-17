// Builds docs/index.html — the illustrated guide (formulas as native MathML, code excerpts, interactive program map).
// Usage:  npm install katex   (anywhere on NODE_PATH)   &&   node docs/guide/build_guide.js
// The output is a single self-contained HTML file: no CDN, no external scripts or fonts.
const fs = require("fs");
const path = require("path");
const katex = require("katex");

const T = String.raw;
const eq = (tex) => '<div class="eq">' + katex.renderToString(tex, { output: "mathml", displayMode: true, throwOnError: true }) + "</div>";
const im = (tex) => katex.renderToString(tex, { output: "mathml", displayMode: false, throwOnError: true });
const esc = (s) => s.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");

// tiny build-time highlighter (C++ / JSON / shell)
const KW = /\b(const|double|int|auto|for|if|else|return|using|static|void|bool|struct|class|std|Eigen|true|false|while|continue|break|switch|case|namespace|inline|constexpr|throw)\b/g;
function hl(src, lang) {
  return src.replace(/\n+$/, "").split("\n").map((line) => {
    let code = line, com = "";
    const m = lang === "cpp" ? line.match(/^(.*?)(\/\/.*)$/) : lang === "bash" ? line.match(/^(.*?)(#.*)$/) : null;
    if (m) { code = m[1]; com = m[2]; }
    let h = esc(code);
    if (lang === "cpp") h = h.replace(KW, "\u0001$1\u0002").replace(/\b(\d+\.?\d*(e[+-]?\d+)?)\b/g, '<span class="n">$1</span>')
      .replace(/\u0001/g, '<span class="k">').replace(/\u0002/g, "</span>");
    if (lang === "json") h = h.replace(/(&quot;|")([^"]*)(&quot;|")(\s*:)/g, '<span class="k">"$2"</span>$4').replace(/:\s*(-?\d+\.?\d*)/g, ': <span class="n">$1</span>');
    return h + (com ? '<span class="c">' + esc(com) + "</span>" : "");
  }).join("\n");
}
const code = (file, src, lang = "cpp") =>
  '<figure class="code"><figcaption>' + esc(file) + "</figcaption><pre><code>" + hl(src, lang) + "</code></pre></figure>";

// ------------------------------------------------------------------ program map
const MAP = [
  { stage: "Inputs", items: [
    ["JSON configuration", "done", "Strict parser: unknown keys are errors; paths relative to the config.", "A reproducible description of one calculation.", "config"],
    ["Aggregate geometry", "done", "15 parametric builders (chains, rings, tubes, sheets, shell on a particle, seamless brickwork tube).", "Positions, transition-dipole directions and molecular frames.", "geometry"],
    ["Optical constants of metals", "done", "refractiveindex.info YAML tables (PCHIP interpolation), Drude model, Kreibig size damping.", "ε(ω) of the particle.", "sphere"],
    ["Transition charges from a geometry", "done", "Hückel π-electron model: q ∝ c(HOMO)·c(LUMO).", "A TrEsp charge file for couplings and dispersion shifts.", "tresp"],
    ["Transition charges from TD-DFT", "partial", "PySCF helper script (CAM-B3LYP, TDA, ESP fit). Works, but is a per-molecule script, not yet a general pipeline.", "Higher-quality TrEsp charges.", "tresp"],
    ["Experimental spectra", "done", "Two- or three-column tables in nm, cm⁻¹ or eV.", "Data for calibration and comparison.", "calibration"],
    ["Atomistic structures (MD)", "planned", "Aggregates on particle surfaces from molecular dynamics.", "Realistic structural disorder.", null],
  ]},
  { stage: "Exciton model", items: [
    ["Frenkel Hamiltonian", "done", "One-exciton block, dense diagonalisation (Eigen).", "Exciton energies, wave functions, transition dipoles, participation factors.", "hamiltonian"],
    ["Excitonic couplings", "done", "Point dipole, extended dipole, atomic transition charges (TrEsp), screening factor.", "The coupling matrix J.", "coupling"],
    ["Disorder", "done", "Gaussian site disorder; spatial correlations via Cholesky factor; orientation/position jitter.", "An ensemble of realisations with error bars.", "disorder"],
    ["Dispersive site shifts", "done", "ΔE = Q·ΣV², calibrated to the Bayliss continuum shift (Megow et al.).", "Environment-dependent site energies — the main source of the J-band red shift in C8S3.", "dispersion"],
    ["Periodic boundary conditions", "planned", "Infinite tubes without end effects.", "Cleaner band widths.", null],
    ["Charge-transfer states", "planned", "Beyond the Frenkel picture.", "Tightly packed dyes.", null],
  ]},
  { stage: "Line shapes and dynamics", items: [
    ["Line-shape function g(t)", "done", "Analytic Matsubara series for the Drude bath; FFT for vibrational modes.", "Input of all line shapes and of modified Redfield theory.", "bath"],
    ["Cumulant spectra", "done", "Second-order cumulant, secular; exchange narrowing through the participation factor.", "Absorption and fluorescence line shapes.", "lineshape"],
    ["Redfield relaxation", "done", "Secular rates between exciton states; detailed balance.", "Kinetics, steady state, lifetime broadening.", "redfield"],
    ["Modified Redfield", "done", "Diagonal exciton–phonon coupling non-perturbative. Practical for N ≲ 300.", "Better rates when λ ~ J.", "modredfield"],
    ["HEOM reference", "done", "Hierarchical equations of motion, Drude bath, small aggregates.", "Numerically exact benchmark for spectra and populations.", "heom"],
    ["Vibronic absorption (Holstein)", "done", "One-particle approximation, one mode per molecule.", "Vibronic band structure of J- and H-aggregates.", "holstein"],
    ["Vibronic fluorescence", "planned", "Two-particle states.", "0–0 / 0–1 emission ratios.", null],
    ["2D and pump–probe spectra", "planned", "Third-order response.", "Comparison with ultrafast experiments.", null],
    ["Raman and Brillouin spectra", "planned", "Vibrational and acoustic response of particles.", "Phonon spectroscopy.", null],
  ]},
  { stage: "Nanoparticle", items: [
    ["Mie theory", "done", "Homogeneous and coated spheres (stable recurrences).", "Exact cross sections of the bare particle.", "sphere"],
    ["Quasi-static multipoles", "done", "Sphere and prolate/oblate spheroids up to l = 120.", "The induced field at every molecule.", "sphere"],
    ["Ligand shell and Feibelman d⊥", "partial", "Layered-sphere boundary conditions. Spheres only; not yet for spheroids.", "Spacer layers and quantum surface response.", "layered"],
    ["Retardation corrections", "done", "Exact Mie factors for l ≤ 8 (sphere); MLWA for spheroids.", "Radiation damping and dynamic depolarisation.", "mlwa"],
    ["Size distribution", "done", "Gauss–Hermite average over the core radius.", "Inhomogeneous plasmon broadening.", "sphere"],
    ["Large particles (R ≳ 20 nm)", "planned", "Radiation of higher multipoles, retardation inside the particle.", "Accurate decay rates for big particles.", null],
    ["Arbitrary shapes, dimers, substrates", "planned", "BEM / DDA back end.", "Gaps and hot spots.", null],
  ]},
  { stage: "Coupled system", items: [
    ["Aggregate + particle response", "done", "Linear response in the exciton basis with the cumulant resolvent.", "Extinction, absorption, scattering of the hybrid.", "response"],
    ["Fast solver", "done", "Woodbury identity for low multipoles, diagonal treatment of high ones.", "~10³ molecules in minutes.", "woodbury"],
    ["Coupling strength g", "done", "From the pole of the dipolar plasmon mode — no peak fitting.", "g, cooperativity, strong-coupling criteria.", "strong"],
    ["Coupled-oscillator analysis", "done", "Branch fit with confidence interval and a 'resolved' flag.", "Anticrossing maps from parameter scans.", "strong"],
    ["Decay near the particle", "done", "Radiation through the particle + absorption in the metal. Accurate for R ≤ 10 nm.", "Purcell factor, quenching, quantum yield, hybrid PL.", "decay"],
  ]},
  { stage: "Calibration and validation", items: [
    ["Monomer calibration", "done", "Joint fit of absorption and PL; moments start, multi-start, Levenberg–Marquardt.", "E₀, bath, static width, μ, radiative lifetime.", "calibration"],
    ["Comparison with experiment", "done", "Scale / baseline / shift fit, residuals, χ², Durbin–Watson.", "A check of what the model misses.", "calibration"],
    ["Benchmarks", "done", "Exact Mie theory for a dipole near a sphere; datasheet extinction of gold.", "Measured accuracy of each block.", "decay"],
    ["Parameter scans", "done", "Any JSON key; maps and observables vs the scanned value.", "Trends and sensitivity.", "config"],
    ["J-band prediction (C8S3)", "partial", "No fit to the aggregate: +100…+200 cm⁻¹ from experiment; band width depends on disorder.", "A test of the whole chain.", "dispersion"],
    ["Uncertainty budget", "planned", "Propagation of input uncertainties.", "Error bars on predictions.", null],
    ["Bayesian calibration", "planned", "Posterior over model parameters, latent variables.", "Systematic model improvement.", null],
  ]},
];

// ------------------------------------------------------------------ methods
const METHODS = [
  { id: "hamiltonian", title: "Frenkel exciton Hamiltonian", status: "done",
    lead: "Each molecule is a two-level system. In the one-exciton manifold the Hamiltonian is a real symmetric N×N matrix; its eigenvectors give the exciton states.",
    eqs: [T`H=\sum_n\left(E_0+\delta_n+\Delta_n\right)|n\rangle\langle n|+\sum_{n\neq m}J_{nm}\,|n\rangle\langle m|`,
      T`\boldsymbol\mu_k=\mu\sum_i c_{ik}\,\hat{\boldsymbol\mu}_i,\qquad s_k=\sum_i c_{ik}^4`],
    notes: "<b>E₀</b> — monomer energy, <b>δ</b> — static disorder, <b>Δ</b> — dispersive shift, <b>J</b> — coupling. The participation factor <i>s</i> (inverse delocalisation length) controls exchange narrowing.",
    code: ["excitons/src/hamiltonian.cpp", `ExcitonStates diagonalize(const Eigen::MatrixXd& H, const Aggregate& agg, double mu_D) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H);
  ExcitonStates s;
  s.energies = es.eigenvalues();
  s.c = es.eigenvectors();                               // c(i, k): site i in state k
  Eigen::MatrixX3d site_mu(N, 3);
  for (Eigen::Index i = 0; i < N; ++i) site_mu.row(i) = mu_D * agg.dipoles[i].transpose();
  s.mu = s.c.transpose() * site_mu;                      // exciton transition dipoles
  s.mu2 = s.mu.rowwise().squaredNorm();
  s.ipr = s.c.array().pow(4).colwise().sum().transpose(); // participation factor s_k
  return s;
}`] },
  { id: "coupling", title: "Excitonic couplings", status: "done",
    lead: "Three levels of detail. The point dipole is valid when molecules are far apart; transition charges keep the shape of the transition density and stay accurate at contact distances.",
    eqs: [T`J^{\mathrm{dip}}_{nm}=f_s\,C\,\mu^2\,\frac{\hat\mu_n\!\cdot\!\hat\mu_m-3(\hat\mu_n\!\cdot\!\hat r)(\hat\mu_m\!\cdot\!\hat r)}{r^3}`,
      T`C=\frac{1}{4\pi\varepsilon_0 hc}=5.034\ \mathrm{cm^{-1}\,nm^{3}\,D^{-2}}`,
      T`J^{\mathrm{TrEsp}}_{nm}=f_s\,C\sum_{a\in n}\sum_{b\in m}\frac{q_a\,q_b}{|\mathbf r_a-\mathbf r_b|}`],
    notes: "The extended dipole is TrEsp with two charges ±μ/l. Charges are rescaled so that |Σqr| equals the measured transition dipole. Try the explorer below.",
    code: ["excitons/src/hamiltonian.cpp", `case CouplingModel::PointDipole: {
  const Vec3 n = rij / r;
  const double kappa = mui.dot(muj) - 3.0 * mui.dot(n) * muj.dot(n);
  return p.screening * units::dd_cm_per_D2_nm3 * p.mu_D * p.mu_D * kappa / (r * r * r);
}
// transition charges: sum over all pairs of atomic charges of the two molecules
for (std::size_t a = 0; a < p.charges.size(); ++a)
  for (std::size_t b = 0; b < p.charges.size(); ++b)
    s += p.charges[a](3) * p.charges[b](3) / (pos[i][a] - pos[j][b]).norm();
J(i, j) = J(j, i) = p.screening * units::dd_cm_per_D2_nm3 * scale * scale * s;`], widget: "w-coupling" },
  { id: "tresp", title: "Transition charges (TrEsp)", status: "partial",
    lead: "Charges come from a file. They can be generated from a molecular geometry with a built-in Hückel calculation, or with the PySCF helper script (TD-DFT, electrostatic-potential fit of the transition density).",
    eqs: [T`q_i\propto c_i^{\mathrm{HOMO}}\,c_i^{\mathrm{LUMO}},\qquad \sum_i q_i=0,\qquad \Big|\sum_i q_i\mathbf r_i\Big|=\mu`],
    notes: "Status: the Hückel generator is part of the program; the TD-DFT route is a script written for one molecule and still has to become a general pipeline.",
    code: ["excitons/src/huckel.cpp", `// Hückel matrix in units of |beta|: heteroatom parameters on the diagonal, bonds off-diagonal
H(k, k) = N ? -(pyrrole ? opt.h_N_pyrrole : opt.h_N_pyridine) : 0.0;
kk = std::exp(-((atoms[i].r - atoms[j].r).norm() - 0.140) / 0.0569);   // C–C: beta depends on bond length
H(k, it->second) = -kk;
...
const Eigen::VectorXd cH = es.eigenvectors().col(R.homo), cL = es.eigenvectors().col(R.homo + 1);
q(k) = std::sqrt(2.0) * cH(k) * cL(k);                                 // HOMO→LUMO transition charge`] },
  { id: "dispersion", title: "Dispersive (van der Waals) site-energy shifts", status: "done",
    lead: "In a dense aggregate each molecule is red-shifted by the polarisable molecules around it. The shift depends on the local packing, so the inner and outer walls of a nanotube absorb at different energies.",
    eqs: [T`\Delta_m=Q\sum_{k\neq m}V_{mk}^2,\qquad V_{mk}=\sum_{a\in m}\sum_{b\in k}\frac{q_a q_b}{\sqrt{r_{ab}^2+a_0^2}}`,
      T`Q\,\Big\langle\sum_k V_{mk}^2\Big\rangle_{\text{random orientations}}=-F\,\frac{n^2-1}{2n^2+1}`,
      T`E_m=E_0^{\text{solvent}}+F\,\frac{n_s^2-1}{2n_s^2+1}+\Delta_m`],
    notes: "After Megow <i>et al.</i>, PCCP 17, 6741 (2015). The softening a₀ avoids a divergence for nearly coincident charges of close neighbours. For C8S3 this term supplies almost the whole J-band shift.",
    code: ["excitons/src/hamiltonian.cpp", `// calibration: random molecular orientations at the same positions give the continuum shift
for (int k = 0; k < K; ++k) {
  for (auto& F : fr) F = Eigen::Quaterniond(g(rng), g(rng), g(rng), g(rng)).normalized().toRotationMatrix();
  avg += sum_v2(agg, place(fr), q, p.cutoff_nm, p.softening_nm).mean();
}
avg /= K;
R.Q_cm_nm2 = R.continuum_cm / avg;
for (std::size_t i = 0; i < agg.size(); ++i) fr[i] = agg.frame(i);            // the real orientations
const Eigen::VectorXd S = sum_v2(agg, place(fr), q, p.cutoff_nm, p.softening_nm);
R.dispersive_cm = R.Q_cm_nm2 * S;
R.shift_cm = R.dispersive_cm.array() + R.to_gas_cm + p.extra_shift_cm;`], widget: "w-shift" },
  { id: "disorder", title: "Disorder and ensemble statistics", status: "done",
    lead: "Spectra are averaged over realisations of the disorder until the standard error is below a tolerance. Observables (peak, width, Stokes shift) are taken from the averaged spectrum; their errors come from a bootstrap.",
    eqs: [T`\langle\delta_i\delta_j\rangle=\sigma^2\,e^{-r_{ij}/\xi},\qquad \boldsymbol\delta=L\,\mathbf z,\quad LL^{\mathsf T}=\mathrm{Cov}`],
    notes: "Structural disorder (random tilts and displacements) rebuilds the coupling matrix in every realisation.",
    code: ["excitons/src/ensemble.cpp", `C(i, j) = model.disorder_sigma_cm * model.disorder_sigma_cm *
          std::exp(-(agg.positions[i] - agg.positions[j]).norm() / model.disorder_corr_length_nm);
Eigen::LLT<Eigen::MatrixXd> llt(C);
chol = llt.matrixL();
...
R.H.diagonal() = Eigen::VectorXd::Constant(z.size(), m.monomer_energy_cm) + (*chol) * z;`] },
  { id: "bath", title: "Bath and line-shape function", status: "done",
    lead: "Each molecule is coupled to its own harmonic bath: an overdamped (Drude) part for the solvent plus underdamped Brownian oscillators for intramolecular modes.",
    eqs: [T`J(\omega)=\frac{2\lambda\gamma\,\omega}{\omega^2+\gamma^2}+\sum_k\frac{2S_k\omega_k^3\gamma_k\,\omega}{(\omega_k^2-\omega^2)^2+\gamma_k^2\omega^2}`,
      T`g(t)=\frac1\pi\int_0^\infty\frac{J(\omega)}{\omega^2}\Big[\coth\frac{\beta\omega}{2}\,(1-\cos\omega t)+i\,(\sin\omega t-\omega t)\Big]\,d\omega`],
    notes: "The Drude part of g(t) is summed analytically (Matsubara series with a tail correction); the mode part is integrated by FFT.",
    code: ["excitons/src/bath.cpp", `double Bath::J_over_w(double w) const {
  double v = 0.0;
  if (lambda_cm > 0.0) v += 2.0 * lambda_cm * wc / (w * w + wc * wc);          // Drude
  for (const auto& m : modes) {                                                 // Brownian oscillators
    const double lk = m.huang_rhys * m.omega_cm;
    const double d = m.omega_cm * m.omega_cm - w * w;
    v += 2.0 * lk * m.omega_cm * m.omega_cm * m.gamma_cm / (d * d + m.gamma_cm * m.gamma_cm * w * w);
  }
  return v;
}`] },
  { id: "lineshape", title: "Cumulant absorption and fluorescence", status: "done",
    lead: "Every exciton state contributes a line whose shape follows from g(t), narrowed by delocalisation and broadened by its population decay.",
    eqs: [T`L^{\mathrm{abs}}_n(\omega)=\frac1\pi\,\mathrm{Re}\int_0^\infty e^{\,i(\omega-\varepsilon_n)t-s_n g(t)-\Gamma_n t/2}\,dt`,
      T`L^{\mathrm{em}}_n(\omega)=\frac1\pi\,\mathrm{Re}\int_0^\infty e^{\,i(\omega-\varepsilon_n+2s_n\lambda)t-s_n g^*(t)-\Gamma_n t/2}\,dt`,
      T`\sigma_{\mathrm{abs}}(\omega)=\frac{\pi\,\omega}{3\varepsilon_0\hbar c\,n}\sum_n|\mu_n|^2L^{\mathrm{abs}}_n(\omega)`],
    notes: "One FFT per state gives the whole line. The same function, kept complex, is the exciton resolvent used in the aggregate + particle response.",
    code: ["excitons/src/spectra.cpp", `const double shift = emission ? 2.0 * s * lambda_ : 0.0;
for (std::size_t j = 0; j < n_; ++j) {
  const double t = dt_ * static_cast<double>(j);
  const std::complex<double> gj = emission ? std::conj(g_[j]) : g_[j];
  std::complex<double> F = std::exp(-s * gj - 0.5 * Gamma_cm * t + std::complex<double>(0.0, shift * t));
  if (j == 0) F *= 0.5;                                  // trapezoid rule
  X[j][0] = F.real();  X[j][1] = F.imag();
}
fftw_execute_dft(static_cast<fftw_plan>(plan_), X, X);   // L(Δ) = dt · Re X / π`] },
  { id: "redfield", title: "Redfield relaxation and kinetics", status: "done",
    lead: "Populations of exciton states obey a Pauli master equation with secular Redfield rates. Detailed balance makes the stationary state Boltzmann.",
    eqs: [T`k_{n\to m}=\sum_i c_{in}^2c_{im}^2\;\tilde C(\varepsilon_n-\varepsilon_m),\qquad \tilde C(\omega)=\frac{2J(|\omega|)\,\mathrm{sgn}\,\omega}{1-e^{-\beta\omega}}`,
      T`\dot P_n=\sum_m\left(k_{m\to n}P_m-k_{n\to m}P_n\right)-(\Gamma^{\mathrm{rad}}_n+\Gamma^{\mathrm{nr}})P_n`],
    notes: "Steady-state populations under pumping weight the fluorescence lines; the quantum yield is ΣΓ<sub>rad</sub>P.",
    code: ["excitons/src/relaxation.cpp", `const Eigen::MatrixXd p = st.c.array().square().matrix();   // p(i,n) = c_in²
const Eigen::MatrixXd overlap = p.transpose() * p;          // Σ_i c_in² c_im²
for (Eigen::Index n = 0; n < N; ++n)
  for (Eigen::Index m = 0; m < N; ++m)
    if (m != n) K(m, n) = overlap(m, n) * bath.rate_spectrum(st.energies(n) - st.energies(m));
for (Eigen::Index n = 0; n < N; ++n) K(n, n) = -K.col(n).sum();`] },
  { id: "modredfield", title: "Modified Redfield theory", status: "done",
    lead: "The diagonal part of the exciton–phonon coupling is treated exactly through g(t); only the off-diagonal part is perturbative. The rate is an overlap of donor fluorescence and acceptor absorption.",
    eqs: [T`k_{M\to N}=2\,\mathrm{Re}\int_0^\infty e^{\,i\Omega t-\Lambda g(t)}\Big\{p\,\ddot g(t)-\big[(q_1-q_2)\,\dot g(t)+2i\lambda q_1\big]^2\Big\}dt`,
      T`\Omega=\varepsilon_M-\varepsilon_N-2\lambda(a-p),\quad \Lambda=a+b-2p,\quad p=\sum_i c_{iM}^2c_{iN}^2,\quad q_1=\sum_i c_{iN}c_{iM}^3`],
    notes: "Checked against detailed balance (0.5 %), the Redfield limit (3 %) and HEOM for a dimer with λ ≈ J.",
    code: ["excitons/src/modified_redfield.cpp", `const double Om = st.energies(M) - st.energies(N) - 2.0 * lam * (a - p);
const double Lam = a + b - 2.0 * p;
for (std::size_t j = 0; j < nt; ++j) {
  const cd e = std::exp(cd(0.0, Om * t) - Lam * g[j]);
  const cd br = (q1 - q2) * g1[j] + cd(0.0, 2.0 * lam * q1);      // g1 = dg/dt, g2 = d²g/dt²
  const double f = (e * (p * g2[j] - br * br)).real();
  acc += (j == 0 ? 0.5 : 1.0) * f;
}
K(N, M) = std::max(0.0, 2.0 * acc * dt);`] },
  { id: "heom", title: "HEOM reference", status: "done",
    lead: "The hierarchical equations of motion solve the same model without the cumulant and Markov approximations. The program uses them as a built-in reference for small aggregates (Drude bath).",
    eqs: [T`\dot\rho_{\mathbf n}=-i[H,\rho_{\mathbf n}]-\sum_{j,k}n_{jk}\nu_k\,\rho_{\mathbf n}-i\sum_{j,k}\big[V_j,\rho_{\mathbf n^{+}_{jk}}\big]-i\sum_{j,k}n_{jk}\big(c_k V_j\rho_{\mathbf n^{-}_{jk}}-c_k^*\rho_{\mathbf n^{-}_{jk}}V_j\big)`],
    notes: "Depth and number of Matsubara terms are parameters; the GUI shows HEOM and the fast theory side by side with the difference." },
  { id: "holstein", title: "Vibronic absorption: Holstein model", status: "done",
    lead: "One intramolecular mode per molecule, one-particle approximation: the excited molecule carries ν vibrational quanta, all others stay in their vibrational ground state.",
    eqs: [T`\langle n\nu|H|n\nu\rangle=E_n-S\omega_0+\nu\omega_0,\qquad \langle n\nu|H|m\nu'\rangle=J_{nm}F_\nu F_{\nu'},\qquad F_\nu=e^{-S/2}\frac{S^{\nu/2}}{\sqrt{\nu!}}`],
    notes: "Reproduces the Poisson progression of the monomer and the growth (J) or loss (H) of the 0–0 band. Vibronic fluorescence needs two-particle states — planned.",
    code: ["excitons/src/vibronic.cpp", `for (Eigen::Index v = 0; v < V; ++v)                            // Franck–Condon factors
  F(v) = std::exp(-0.5 * p.huang_rhys + 0.5 * v * std::log(std::max(p.huang_rhys, 1e-300)) - 0.5 * std::lgamma(v + 1.0));
const Eigen::MatrixXd FF = F * F.transpose();
...
Hv(n * V + v, n * V + v) = H(n, n) - p.huang_rhys * p.omega_cm + double(v) * p.omega_cm;
Hv.block(n * V, m * V, V, V) = H(n, m) * FF;`] },
  { id: "sphere", title: "Nanoparticle: multipole response", status: "done",
    lead: "The field that the particle adds at each molecule is expanded in multipoles. Close to the surface many orders are needed (up to l = 120); the bare particle's cross sections come from full Mie theory.",
    eqs: [T`\hat D_l=-\frac{l\,(\varepsilon-\varepsilon_m)}{l\,\varepsilon+(l+1)\,\varepsilon_m},\qquad \alpha=-\hat D_1R^3=R^3\frac{\varepsilon-\varepsilon_m}{\varepsilon+2\varepsilon_m}`,
      T`W_{ij}(\omega)=\frac{C\mu^2}{\varepsilon_m}\,(\hat\mu_i\!\cdot\!\nabla)(\hat\mu_j\!\cdot\!\nabla')\,K_{\mathrm{ind}}(\mathbf r_i,\mathbf r_j;\omega)`],
    notes: "For spheroids the same expansion is written in spheroidal coordinates with Legendre functions of the second kind. A Gauss–Hermite average over the radius models a size distribution.",
    code: ["excitons/src/plexciton.cpp", `// sphere: D_l depends only on l, so W = pref · Σ_l D_l(ω) · G_l with precomputed N×N matrices G_l
for (int l = 1; l <= particle_->l_max(); ++l) {
  const cd Dl = d(t0) / terms[t0].w;
  Wr.noalias() += Dl.real() * Gl_[l - 1];
  Wi.noalias() += Dl.imag() * Gl_[l - 1];
}
F.W = pref * (Wr.cast<cd>() + cd(0.0, 1.0) * Wi.cast<cd>());`] },
  { id: "layered", title: "Ligand shell and Feibelman surface response", status: "partial",
    lead: "A dielectric shell around the metal core and the quantum spill-out of the electron density both enter as boundary conditions of the multipole problem.",
    eqs: [T`\varepsilon_{\mathrm{in}}\,\partial_r\phi_{\mathrm{in}}=\varepsilon_{\mathrm{out}}\,\partial_r\phi_{\mathrm{out}},\qquad \phi_{\mathrm{out}}-\phi_{\mathrm{in}}=-d_\perp\,(\partial_r\phi_{\mathrm{out}}-\partial_r\phi_{\mathrm{in}})`],
    notes: "Radial functions are scaled by the outer radius — otherwise the linear system becomes singular for l ≳ 7. Status: spheres only; spheroids are planned.",
    code: ["plasmonics/src/layered.cpp", `const double rp = std::pow(R / bR, l), rq = std::pow(bR / R, l + 1);   // scaled r^l and r^-(l+1)
const double dp = L / R * rp, dq = -(L + 1) / R * rq;
// (1) eps_in dphi_in - eps_out dphi_out = 0      (2) phi_out - phi_in + d (dphi_out - dphi_in) = 0
A(e1, col) += eps[b] * der;
A(e2, col) += -val - d * der;`] },
  { id: "mlwa", title: "Retardation corrections", status: "done",
    lead: "The quasi-static response is corrected for radiation damping and dynamic depolarisation: exact Mie coefficients for spheres, the modified long-wavelength approximation for spheroids.",
    eqs: [T`f_l=\frac{a_l^{\mathrm{Mie}}(x)}{a_l^{\mathrm{qs}}(x)},\qquad a_l^{\mathrm{qs}}=\frac{i\,(l+1)\,x^{2l+1}}{l\,(2l+1)!!\,(2l-1)!!}\,\hat D_l`,
      T`\alpha'=\frac{\alpha}{1-\dfrac{k^2}{a}\alpha-\dfrac{2i}{3}k^3\alpha}`],
    code: ["excitons/src/plexciton.cpp", `// MLWA (Meier–Wokaun): alpha' = alpha / (1 - (k²/a) alpha - i (2/3) k³ alpha)
a_dyn(a) = a_qs(a) / (1.0 - F.k * F.k / ax(a) * a_qs(a) - cd(0.0, 2.0 / 3.0) * std::pow(F.k, 3) * a_qs(a));`] },
  { id: "response", title: "Aggregate + particle: linear response", status: "done",
    lead: "Molecules and particle are solved together. The exciton block uses the same cumulant resolvent as the bare-aggregate spectrum, so both spectra are consistent by construction.",
    eqs: [T`\alpha_{\mathrm{tot}}(\omega)=\alpha_{\mathrm{NP}}(\omega)+c\,M^{\mathsf T}\big[\mathrm{diag}\,a_n(\omega)+U^{\mathsf T}W(\omega)U\big]^{-1}M,\qquad a_n=-\frac{i}{X_n(\omega)}`,
      T`\sigma_{\mathrm{ext}}=4\pi k\,\mathrm{Im}\,\tfrac13\mathrm{Tr}\,\alpha_{\mathrm{tot}},\qquad \sigma_{\mathrm{sca}}=\tfrac{8\pi}{3}k^4\,\tfrac13\sum_{ab}|\alpha_{ab}|^2,\qquad \sigma_{\mathrm{abs}}=\sigma_{\mathrm{ext}}-\sigma_{\mathrm{sca}}`],
    notes: "X<sub>n</sub>(ω) is the complex cumulant response of exciton n; M contains the local-field factors L<sub>i</sub><sup>T</sup>μ̂<sub>i</sub> projected on the exciton states.",
    code: ["excitons/src/plexciton.cpp", `// in the exciton basis: A = diag(a_n) + Uᵀ W U,  M' = Uᵀ M
Eigen::MatrixXcd A = U.transpose().cast<cd>() * F.W * U.cast<cd>();
for (Eigen::Index n = 0; n < N; ++n) A(n, n) += diag[rr][n * nw + q];
X = A.partialPivLu().solve(Mx);
const Eigen::Matrix3cd a_tot = a_np + c_mol * (Mx.transpose() * X);
// cross(a, ext, sca) for a = a_tot:
ext = 4.0 * pi * F.k * a.trace().imag() / 3.0;
sca = 8.0 * pi / 3.0 * std::pow(F.k, 4) * a.cwiseAbs2().sum() / 3.0;`] },
  { id: "woodbury", title: "Fast solver for large aggregates", status: "done",
    lead: "Low multipoles couple all molecules but have low rank; high multipoles are short-ranged and nearly diagonal in the exciton basis. This turns an N×N solve per frequency into a small one.",
    eqs: [T`A=A_0+B\,\delta\,B^{\mathsf T},\qquad X=Y-Z\,(I+\delta B^{\mathsf T}Z)^{-1}\,\delta B^{\mathsf T}Y,\qquad Y=A_0^{-1}M,\ \ Z=A_0^{-1}B`],
    notes: "Error against the exact solver for 300 molecules: 0.7 % (l ≤ 8 exact) and 0.16 % (l ≤ 16).",
    code: ["excitons/src/plexciton.cpp", `const Eigen::MatrixXcd Y = inv.asDiagonal() * Mx, Z = inv.asDiagonal() * B;
Eigen::MatrixXcd S = dl.asDiagonal() * (B.transpose() * Z);
S.diagonal().array() += 1.0;
X = Y - Z * S.partialPivLu().solve(dl.asDiagonal() * (B.transpose() * Y));`] },
  { id: "decay", title: "Decay near the particle: Purcell factor and quantum yield", status: "done",
    lead: "A molecule next to the metal radiates together with the dipole it induces in the particle, and loses energy to absorption in the metal. Their balance decides between enhancement and quenching.",
    eqs: [T`\Gamma_{\mathrm{tot}}=\Gamma_0\,\big|L^{\mathsf T}\hat\mu\big|^2+\Gamma_{\mathrm{abs}},\qquad \mathrm{QY}=\frac{\Gamma_0\,|L^{\mathsf T}\hat\mu|^2}{\Gamma_{\mathrm{tot}}+\Gamma_{\mathrm{nr}}}`,
      T`\Gamma_0=\frac{n\,\omega^3\mu^2}{3\pi\varepsilon_0\hbar c^3}`],
    notes: "Validated against the exact Mie solution for a dipole near a sphere: 1–3 % for R ≤ 10 nm; the error grows for R ≳ 20 nm (see Validation).",
    code: ["excitons/src/plexciton.cpp", `const Eigen::Vector3cd me = F.L[i].transpose() * agg_.dipoles[i].cast<cd>();
o.rad_enhancement = me.squaredNorm();                              // |Lᵀ mu|²: emission through the particle
o.nonrad_to_np_ps = std::max(0.0, (-2.0 * Wa(ii).imag() + 2.0 * pref * ret) * units::rad_ps_per_cm);
const double grad = gamma0_ps * o.rad_enhancement;
o.quantum_yield = grad / (grad + o.nonrad_to_np_ps + gamma_nr_ps);`] },
  { id: "strong", title: "Coupling strength and strong-coupling criteria", status: "done",
    lead: "The collective coupling g between the molecular layer and the dipolar plasmon is computed from the plasmon pole and the mode field at the molecules — not from fitted peak positions, which overestimate g when the branches are not resolved.",
    eqs: [T`\alpha_a(\omega)\approx\frac{A}{\omega_p^2-\omega^2-i\gamma_p\omega}+B,\qquad g^2=\frac{A}{2\omega_p}\,\mathrm{Tr}\,K_a`,
      T`\omega_\pm=\frac{\omega_p+\omega_x}{2}-\frac{i(\gamma_p+\gamma_x)}{4}\pm\sqrt{g^2+\Big[\frac{\omega_p-\omega_x}{2}-\frac{i(\gamma_p-\gamma_x)}{4}\Big]^2}`,
      T`2g>\frac{\gamma_p+\gamma_x}{2}\ \ (\text{visible splitting}),\qquad g>\frac{|\gamma_p-\gamma_x|}{4}\ \ (\text{anticrossing}),\qquad C=\frac{4g^2}{\gamma_p\gamma_x}`],
    code: ["excitons/app/scan.cpp", `void coupled_branches(double wp, double gp, double wx, double gx, double g, double& lower, double& upper) {
  const cd mean = 0.5 * (wp + wx) - cd(0.0, 0.25 * (gp + gx));
  const cd half = 0.5 * (wp - wx) - cd(0.0, 0.25 * (gp - gx));
  const cd root = std::sqrt(g * g + half * half);
  const double a = (mean + root).real(), b = (mean - root).real();
  lower = std::min(a, b);
  upper = std::max(a, b);
}`], widget: "w-anticross" },
  { id: "calibration", title: "Calibration on monomer spectra", status: "done",
    lead: "Before an aggregate is calculated, the molecule itself is calibrated: absorption and fluorescence of the monomer are fitted together, which fixes the transition energy, the bath and the static width. The transition dipole follows from the absolute extinction.",
    eqs: [T`A(\omega)\propto\omega\,\big[L^{\mathrm{abs}}\otimes G_\sigma\big](\omega),\qquad F(\omega)\propto\omega^3\,\big[L^{\mathrm{em}}\otimes G_\sigma\big](\omega)`,
      T`\mu^2=\frac{3\varepsilon_0\hbar c\,n}{\pi}\int\frac{\sigma(\tilde\nu)}{\tilde\nu}\,d\tilde\nu,\qquad k_r=2.880\times10^{-9}\,n^2\,\langle\tilde\nu^{-3}\rangle^{-1}\!\int\frac{\epsilon(\tilde\nu)}{\tilde\nu}\,d\tilde\nu`],
    notes: "Static broadening is applied in the time domain (a factor exp(−σ²t²/2)), so the model is evaluated directly at the data points. Without a fluorescence spectrum the bath is not determined — always supply one.",
    code: ["excitons/src/monomer_fit.cpp", `const double env = std::exp(-0.5 * sg * sg * t * t);                       // static Gaussian broadening
fa[j] = wj * env * std::exp(-g[j]);                                         // absorption kernel
fe[j] = wj * env * std::exp(std::complex<double>(0.0, 2.0 * lam * t) - std::conj(g[j]));   // emission kernel
...
abs[i] = x * sa.real();             // A(w) ∝ w · L_abs
em[i]  = x * x * x * se.real();     // F(w) ∝ w³ · L_em`] },
];
const EXTRA = { geometry: "Geometry builders live in excitons/src/geometry.cpp; run ssp_aggregate --list-geometries.", config: "See the Configuration section." };

// ------------------------------------------------------------------ static content
const STATUS = { done: ["Implemented", "✓"], partial: ["Partly implemented", "◐"], planned: ["Not implemented yet", "○"] };
const badge = (s) => '<span class="badge ' + s + '"><i>' + STATUS[s][1] + "</i>" + STATUS[s][0] + "</span>";

const methodsHtml = METHODS.map((m) => `
<article class="method" id="m-${m.id}">
  <header><h3>${m.title}</h3>${badge(m.status)}</header>
  <p class="lead">${m.lead}</p>
  <div class="cols">
    <div class="math">${m.eqs.map(eq).join("")}${m.notes ? '<p class="note">' + m.notes + "</p>" : ""}</div>
    ${m.code ? '<div class="src">' + code(m.code[0], m.code[1]) + "</div>" : ""}
  </div>
  ${m.widget ? '<div class="widget-slot" data-widget="' + m.widget + '"></div>' : ""}
</article>`).join("\n");

const cfg = (title, rows) => `<details><summary>${title}</summary><table><thead><tr><th>Key</th><th>Default</th><th>Meaning</th></tr></thead><tbody>${
  rows.map((r) => `<tr><td><code>${r[0]}</code></td><td>${r[1]}</td><td>${r[2]}</td></tr>`).join("")}</tbody></table></details>`;
const configHtml = [
  cfg("<code>configurations[]</code> — which aggregates to calculate", [
    ["name", "—", "label of the shape; also the name of its results folder"],
    ["geometry_type", "—", "one of the 15 builders (<code>ssp_aggregate --list-geometries</code>)"],
    ["geometry_params", "builder defaults", "numbers specific to the builder (sizes in nm, angles in degrees)"],
    ["placement.shift_nm, euler_deg, center", "0, 0, false", "move and rotate the aggregate relative to the particle at the origin"]]),
  cfg("<code>exciton</code> — the molecule, couplings and disorder", [
    ["monomer_energy_cm", "17000", "vertical transition energy of the isolated molecule, cm⁻¹"],
    ["transition_dipole_D", "10", "transition dipole, debye"],
    ["disorder_sigma_cm", "60", "standard deviation of the site energies"],
    ["disorder_corr_length_nm", "0", "correlation length ξ of the site energies (0 — uncorrelated)"],
    ["orientation_jitter_deg, position_jitter_nm", "0", "structural disorder, drawn anew in every realisation"],
    ["relaxation", "redfield", "<code>redfield</code> or <code>modified_redfield</code>"],
    ["coupling.model", "point_dipole", "<code>point_dipole</code>, <code>extended_dipole</code>, <code>transition_charges</code>, <code>legacy_J0</code>"],
    ["coupling.screening", "1/n²", "factor f<sub>s</sub> multiplying all couplings"],
    ["coupling.charges_file", "—", "TrEsp charges: columns x y z q (nm, molecular frame, x along the dipole)"],
    ["dispersion.enabled", "false", "switch on environment-dependent site shifts"],
    ["dispersion.F_eV, n_aggregate, n_solvent", "1.43, 1.78, n_medium", "Bayliss factor, refractive index of the dye layer and of the monomer's solvent"],
    ["dispersion.softening_nm, cutoff_nm", "0.3, 3.0", "charge softening a₀ and summation radius"],
    ["dispersion.use_coupling_charges", "false", "use the TrEsp charges of the coupling for V<sub>mk</sub>"]]),
  cfg("<code>bath</code>, <code>decay</code>, <code>pump</code> — line shapes and emission", [
    ["bath.type", "drude_lorentz", "or <code>ohmic_exponential</code>"],
    ["bath.lambda_cm, cutoff_cm", "50, 100", "reorganisation energy and cut-off of the overdamped part"],
    ["bath.temperature_K", "300", "temperature"],
    ["bath.modes[]", "[]", "vibrational modes: <code>omega_cm</code>, <code>huang_rhys</code>, <code>gamma_cm</code>"],
    ["decay.monomer_lifetime_ns", "0", "radiative lifetime of the monomer (0 — from μ and ω)"],
    ["decay.nonradiative_ps", "0", "intrinsic non-radiative rate, ps⁻¹"],
    ["pump.type, excitation_cm", "oscillator_strength", "how the states are populated: by oscillator strength or by a laser line"]]),
  cfg("<code>spectrum</code>, <code>ensemble</code>, <code>dynamics</code> — numerics", [
    ["spectrum.wmin_cm, wmax_cm, nw", "15500, 18000, 1800", "spectral window and number of points"],
    ["spectrum.t_max_ps", "100", "time window of the line-shape integrals"],
    ["ensemble.min_realizations, max_realizations", "20, 400", "limits of the disorder average"],
    ["ensemble.rel_tolerance", "0.01", "stop when the standard error / peak is below this value"],
    ["dynamics.enabled, t_max_ps, n_times, initial", "true, 20, 200, brightest", "population kinetics"]]),
  cfg("<code>nanoparticle</code> — the particle and the hybrid calculation", [
    ["enabled", "false", "include a particle centred at the origin"],
    ["shape", "sphere", "<code>sphere</code>, <code>prolate</code> (rod), <code>oblate</code> (disk)"],
    ["semi_axis_xy_nm, semi_axis_z_nm", "15, 15", "semi-axes of the metal core"],
    ["material_file | drude", "—", "refractiveindex.info YAML file, or Drude parameters <code>eps_inf, omega_p_eV, gamma_eV</code>"],
    ["n_medium", "1.33", "refractive index of the surrounding medium"],
    ["l_max", "40", "highest multipole; the log warns when it is too low"],
    ["shell_thickness_nm, shell_n", "0, 1.45", "ligand / spacer shell (sphere)"],
    ["d_perp_nm, d_perp_im_nm", "0, 0", "Feibelman parameter d⊥ (sphere)"],
    ["kreibig_A", "0", "size-dependent damping, γ = γ<sub>bulk</sub> + A·v<sub>F</sub>/R"],
    ["size_sigma_rel, size_nodes", "0, 5", "relative spread of the core radius and number of quadrature nodes"],
    ["fast_l", "auto", "<code>auto</code>, <code>0</code> (exact) or the highest exactly treated multipole"],
    ["emission", "true", "hybrid PL, quantum yield, Purcell factors"]]),
  cfg("<code>holstein</code>, <code>heom</code>, <code>scan</code>, <code>monomer_fit</code> — extras", [
    ["holstein.enabled, omega_cm, huang_rhys, nu_max", "false, 1300, 0.6, 4", "vibronic absorption"],
    ["heom.enabled, depth, matsubara, max_sites", "false, 4, 1, 10", "exact reference for small aggregates"],
    ["scan.parameter, from, to, steps", "—", "scan any JSON path, e.g. <code>exciton.disorder_sigma_cm</code>"],
    ["monomer_fit.absorption_file, emission_file", "—", "calibration run: spectra, units, start values, fit flags"]]),
].join("\n");

const html = `<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>SSP Nano — illustrated guide</title>
<meta name="description" content="Formulas, code and an interactive map of SSP Nano: exciton spectra of molecular aggregates and plasmon–exciton hybrids.">
<style>
:root{--bg:#fbfbf9;--panel:#fff;--panel2:#f3f3ef;--ink:#1c1c1a;--ink2:#55544e;--line:#deddd6;--accent:#2a78d6;--accent2:#eb6834;
--done:#0f7a52;--done-bg:#dff3ea;--partial:#8a5a00;--partial-bg:#fbecc8;--planned:#5b5a66;--planned-bg:#e8e7ee;--code:#f5f4ef;--shadow:0 1px 2px rgba(0,0,0,.06),0 6px 20px rgba(0,0,0,.05)}
@media (prefers-color-scheme:dark){:root:not([data-theme=light]){--bg:#161615;--panel:#1f1f1d;--panel2:#262624;--ink:#f0efec;--ink2:#b5b4aa;--line:#3a3a36;--accent:#5a9bea;--accent2:#f08a5d;
--done:#5fd3a4;--done-bg:#143528;--partial:#f0c060;--partial-bg:#3b2f10;--planned:#b9b8c6;--planned-bg:#2e2d36;--code:#1b1b19;--shadow:none}}
:root[data-theme=dark]{--bg:#161615;--panel:#1f1f1d;--panel2:#262624;--ink:#f0efec;--ink2:#b5b4aa;--line:#3a3a36;--accent:#5a9bea;--accent2:#f08a5d;
--done:#5fd3a4;--done-bg:#143528;--partial:#f0c060;--partial-bg:#3b2f10;--planned:#b9b8c6;--planned-bg:#2e2d36;--code:#1b1b19;--shadow:none}
*{box-sizing:border-box}html{scroll-behavior:smooth;scroll-padding-top:16px}
body{margin:0;background:var(--bg);color:var(--ink);font:16px/1.6 system-ui,-apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}
a{color:var(--accent)}code{font:13.5px/1.5 ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;background:var(--panel2);padding:.08em .35em;border-radius:4px}
.layout{display:grid;grid-template-columns:230px minmax(0,1fr);max-width:1320px;margin:0 auto}
nav{position:sticky;top:0;align-self:start;height:100vh;overflow:auto;padding:22px 14px 22px 20px;border-right:1px solid var(--line)}
nav .brand{font-weight:700;font-size:18px;margin-bottom:2px}nav .sub{color:var(--ink2);font-size:12.5px;margin-bottom:14px}
nav a{display:block;color:var(--ink2);text-decoration:none;padding:4px 8px;border-radius:6px;font-size:14px}nav a:hover{background:var(--panel2);color:var(--ink)}
nav .grp{margin:12px 0 4px 8px;font-size:11px;letter-spacing:.08em;text-transform:uppercase;color:var(--ink2)}
nav button{margin-top:14px;font:inherit;font-size:13px;padding:6px 10px;border:1px solid var(--line);background:var(--panel);color:var(--ink);border-radius:8px;cursor:pointer}
main{padding:28px 32px 80px;min-width:0}
h1{font-size:38px;line-height:1.15;margin:0 0 10px;letter-spacing:-.02em}h2{font-size:26px;margin:56px 0 8px;letter-spacing:-.01em}h3{font-size:19px;margin:0}
.tag{color:var(--ink2);font-size:18px;max-width:60ch;margin:0 0 18px}
.hero{display:grid;grid-template-columns:minmax(0,1fr) minmax(0,1.25fr);gap:28px;align-items:center}
.hero img,.shots img{width:100%;border:1px solid var(--line);border-radius:10px;box-shadow:var(--shadow);display:block}
.kpis{display:flex;flex-wrap:wrap;gap:10px;margin-top:8px}.kpi{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:8px 14px}.kpi b{display:block;font-size:20px}.kpi span{font-size:12.5px;color:var(--ink2)}
.badge{display:inline-flex;align-items:center;gap:6px;font-size:12px;font-weight:600;padding:2px 9px 2px 7px;border-radius:999px;white-space:nowrap}
.badge i{font-style:normal}.badge.done{color:var(--done);background:var(--done-bg)}.badge.partial{color:var(--partial);background:var(--partial-bg)}.badge.planned{color:var(--planned);background:var(--planned-bg)}
.chips{display:flex;flex-wrap:wrap;gap:8px;margin:14px 0}.chip{font:inherit;font-size:13.5px;padding:5px 12px;border-radius:999px;border:1px solid var(--line);background:var(--panel);color:var(--ink);cursor:pointer}
.chip[aria-pressed=true]{background:var(--ink);color:var(--bg);border-color:var(--ink)}
.map{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:14px}
.stage{background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:12px}
.stage h4{margin:0 0 8px;font-size:13px;letter-spacing:.06em;text-transform:uppercase;color:var(--ink2);display:flex;justify-content:space-between}
.node{display:flex;align-items:center;gap:8px;width:100%;text-align:left;font:inherit;font-size:14px;padding:7px 9px;margin:4px 0;border-radius:8px;border:1px solid var(--line);background:var(--bg);color:var(--ink);cursor:pointer}
.node:hover,.node[aria-expanded=true]{border-color:var(--accent)}.node .dot{flex:none;width:20px;height:20px;border-radius:50%;display:grid;place-items:center;font-size:12px;font-weight:700}
.node.done .dot{background:var(--done-bg);color:var(--done)}.node.partial .dot{background:var(--partial-bg);color:var(--partial)}.node.planned .dot{background:var(--planned-bg);color:var(--planned)}
.node.planned{border-style:dashed;color:var(--ink2)}.node[hidden]{display:none}
.detail{margin-top:14px;background:var(--panel);border:1px solid var(--accent);border-radius:12px;padding:16px 18px;min-height:96px}
.detail h4{margin:0 0 6px;font-size:17px;display:flex;gap:10px;align-items:center;flex-wrap:wrap}.detail dl{display:grid;grid-template-columns:110px 1fr;gap:4px 12px;margin:8px 0 0}.detail dt{color:var(--ink2);font-size:13.5px}.detail dd{margin:0}
.flow{display:flex;flex-wrap:wrap;gap:6px;align-items:center;margin:10px 0 2px;color:var(--ink2);font-size:13.5px}.flow b{background:var(--panel2);border-radius:6px;padding:2px 8px;color:var(--ink);font-weight:600}
.method{background:var(--panel);border:1px solid var(--line);border-radius:14px;padding:18px 20px;margin:18px 0;box-shadow:var(--shadow)}
.method header{display:flex;justify-content:space-between;gap:12px;align-items:center;flex-wrap:wrap}.lead{color:var(--ink2);margin:8px 0 6px;max-width:80ch}
.cols{display:grid;grid-template-columns:minmax(0,1fr);gap:12px}.math{max-width:860px}
.eq{overflow-x:auto;padding:8px 2px;font-size:17px}math{font-family:"STIX Two Math","Cambria Math","Latin Modern Math","Noto Sans Math",serif}
.note{font-size:14px;color:var(--ink2);margin:6px 0 0}
figure.code{margin:0;border:1px solid var(--line);border-radius:10px;overflow:hidden;background:var(--code)}
figure.code figcaption{font:12px ui-monospace,Menlo,Consolas,monospace;color:var(--ink2);padding:6px 12px;border-bottom:1px solid var(--line);background:var(--panel2)}
pre{margin:0;padding:12px;overflow-x:auto;font:12.5px/1.55 ui-monospace,SFMono-Regular,Menlo,Consolas,monospace}pre code{background:none;padding:0;font:inherit}
.k{color:var(--accent)}.n{color:var(--accent2)}.c{color:var(--ink2);font-style:italic}
.widget{margin-top:14px;border-top:1px dashed var(--line);padding-top:14px;display:grid;grid-template-columns:minmax(0,300px) minmax(0,1fr);gap:18px;align-items:center}
.widget h5{margin:0 0 8px;font-size:14px;grid-column:1/-1}.widget label{display:grid;grid-template-columns:1fr auto;font-size:13.5px;color:var(--ink2);margin-top:8px}
.widget label output{color:var(--ink);font-variant-numeric:tabular-nums}.widget input[type=range]{grid-column:1/-1;width:100%;accent-color:var(--accent)}
.widget canvas{width:100%;height:230px;background:var(--bg);border:1px solid var(--line);border-radius:10px}.readout{margin-top:10px;font-size:14px}.readout b{font-variant-numeric:tabular-nums}
.shots{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:18px}.shots figure{margin:0}.shots figcaption{font-size:13.5px;color:var(--ink2);margin-top:6px}
table{border-collapse:collapse;width:100%;font-size:14px;margin:10px 0}th,td{border-bottom:1px solid var(--line);padding:7px 10px;text-align:left;vertical-align:top}th{color:var(--ink2);font-weight:600;font-size:12.5px;text-transform:uppercase;letter-spacing:.04em}
.tablewrap{overflow-x:auto}details{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:8px 14px;margin:8px 0}summary{cursor:pointer;font-weight:600}
.tabs{display:flex;gap:6px;margin:12px 0 0}.tab{font:inherit;font-size:14px;padding:6px 14px;border:1px solid var(--line);border-bottom:none;border-radius:8px 8px 0 0;background:var(--panel2);color:var(--ink2);cursor:pointer}
.tab[aria-selected=true]{background:var(--panel);color:var(--ink);font-weight:600}.tabpanel{border:1px solid var(--line);border-radius:0 10px 10px 10px;background:var(--panel);padding:14px 16px}.tabpanel[hidden]{display:none}
.road{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:12px}.road div{background:var(--panel);border:1px dashed var(--line);border-radius:10px;padding:12px 14px}.road b{display:block}
footer{margin-top:60px;color:var(--ink2);font-size:13.5px;border-top:1px solid var(--line);padding-top:16px}
@media (max-width:980px){.layout{grid-template-columns:1fr}nav{position:static;height:auto;border-right:none;border-bottom:1px solid var(--line);padding:14px 16px}nav .links{display:flex;flex-wrap:wrap;gap:2px}nav .grp{display:none}
main{padding:20px 16px 60px}.hero,.cols,.widget,.shots,.road{grid-template-columns:1fr}.map{grid-template-columns:1fr}h1{font-size:30px}}
@media (min-width:981px) and (max-width:1180px){.map{grid-template-columns:repeat(2,minmax(0,1fr))}}
</style>
</head>
<body>
<div class="layout">
<nav>
  <div class="brand">SSP Nano</div><div class="sub">illustrated guide · v0.7.3</div>
  <div class="links">
  <a href="#top">Overview</a><a href="#map">Program map</a>
  <div class="grp">Methods</div>
  ${METHODS.map((m) => `<a href="#m-${m.id}">${m.title.replace(/:.*/, "")}</a>`).join("")}
  <div class="grp">Using it</div>
  <a href="#gallery">Interface</a><a href="#config">Configuration</a><a href="#build">Build</a><a href="#validation">Validation</a><a href="#roadmap">Roadmap</a>
  </div>
  <button id="theme" type="button" aria-label="Switch colour theme">◐ Theme</button>
</nav>
<main id="top">
<section class="hero">
  <div>
    <h1>Exciton spectra of molecular aggregates and plasmon–exciton hybrids</h1>
    <p class="tag">SSP Nano calculates absorption, fluorescence, exciton kinetics and quantum yields of J- and H-aggregates of any shape — alone or next to a metal nanoparticle — and checks itself against exact theory and experiment.</p>
    <div class="kpis">
      <div class="kpi"><b>15</b><span>aggregate shapes</span></div>
      <div class="kpi"><b>61</b><span>automated tests</span></div>
      <div class="kpi"><b>~10³</b><span>molecules near a particle</span></div>
      <div class="kpi"><b>C++20 · Qt 6</b><span>CLI and GUI</span></div>
    </div>
  </div>
  <img src="screenshots/gui_spectra_vs_experiment.png" alt="Main window of SSP Nano: predicted spectrum of a C8S3 nanotube compared with the measured spectrum, with a residuals panel">
</section>

<h2 id="map">Program map</h2>
<p class="lead">What the program consists of, which method does each job, what it produces — and what is not there yet. Select a block to see the details; use the filter to show only one status.</p>
<div class="flow"><b>structure</b>→<b>Hamiltonian</b>→<b>exciton states</b>→<b>line shapes + rates</b>→<b>(+ particle response)</b>→<b>spectra, kinetics, yields</b>→<b>comparison with experiment</b></div>
<div class="chips" role="group" aria-label="Filter by status">
  <button class="chip" data-f="all" aria-pressed="true">All <span id="c-all"></span></button>
  <button class="chip" data-f="done" aria-pressed="false">✓ Implemented <span id="c-done"></span></button>
  <button class="chip" data-f="partial" aria-pressed="false">◐ Partly <span id="c-partial"></span></button>
  <button class="chip" data-f="planned" aria-pressed="false">○ Not yet <span id="c-planned"></span></button>
</div>
<div class="map" id="mapgrid"></div>
<div class="detail" id="detail" aria-live="polite"><h4>Select a block</h4><p class="note">Each block shows the method, what it gives, and a link to the formulas and code.</p></div>

<h2 id="methods">Methods: formulas and the code behind them</h2>
<p class="lead">Units: energies in cm⁻¹, lengths in nm, dipoles in debye, ħ = 1. Code excerpts are abridged from the files named above each block (comments translated to English); paths are relative to <code>software/</code>.</p>
${methodsHtml}

<h2 id="gallery">The interface</h2>
<p class="lead">Parameters and a 3D preview on the left, results on the right. Every plot can be zoomed, exported (PNG, SVG, CSV) and compared with measured data.</p>
<div class="shots">
  <figure><img loading="lazy" src="screenshots/gui_calibration.png" alt="Calibration tab"><figcaption><b>Calibration.</b> Monomer absorption and fluorescence, fitted model, residuals and parameters with errors.</figcaption></figure>
  <figure><img loading="lazy" src="screenshots/gui_scan_hybrid.png" alt="Scan tab"><figcaption><b>Scan.</b> The exciton energy is tuned across the plasmon of a gold sphere; coupling analysis below the map.</figcaption></figure>
  <figure><img loading="lazy" src="screenshots/gui_aggregate_np.png" alt="Aggregate plus nanoparticle tab"><figcaption><b>Aggregate + particle.</b> Contribution of the molecular layer to extinction, absorption and scattering.</figcaption></figure>
  <figure><img loading="lazy" src="screenshots/gui_near_field.png" alt="Near field tab"><figcaption><b>Near field.</b> Field enhancement, quenching rate and quantum yield for every molecule.</figcaption></figure>
  <figure><img loading="lazy" src="screenshots/gui_heom.png" alt="HEOM reference tab"><figcaption><b>HEOM reference.</b> The fast theory against the numerically exact calculation.</figcaption></figure>
  <figure><img loading="lazy" src="demo/figures/fig2_nanotube_levels.png" alt="Predicted and measured J-band of the C8S3 nanotube"><figcaption><b>Example study.</b> J-band of the C8S3 nanotube predicted at four levels of the model, without fitting the aggregate spectrum.</figcaption></figure>
</div>

<h2 id="config">Configuration</h2>
<p class="lead">One JSON file describes a calculation. Unknown keys are errors; keys that start with an underscore are comments; file paths are relative to the config. The GUI edits the same file.</p>
${code("examples/demo_c8s3/2d_nanotube_tddft.json (abridged)", `{
  "output_dir": "results/2_nanotube",
  "configurations": [ { "name": "tddft_charges", "geometry_type": "brickwork_tube", "geometry_params": {} } ],
  "exciton": {
    "monomer_energy_cm": 19906.5, "transition_dipole_D": 11.4, "disorder_sigma_cm": 124.5,
    "coupling":   { "model": "transition_charges", "charges_file": "data/c8s3_tddft_tresp_S1.txt", "screening": 0.32 },
    "dispersion": { "enabled": true, "F_eV": 1.43, "n_aggregate": 1.78, "n_solvent": 1.329, "use_coupling_charges": true }
  },
  "bath": { "type": "drude_lorentz", "lambda_cm": 764.0, "cutoff_cm": 572.6, "temperature_K": 295,
            "modes": [ { "omega_cm": 1258.6, "huang_rhys": 0.1957, "gamma_cm": 30 } ] },
  "spectrum": { "wmin_cm": 14500, "wmax_cm": 21500, "nw": 1400 }
}`, "json")}
${configHtml}

<h2 id="build">Build and run</h2>
<div class="tabs" role="tablist">
  <button class="tab" role="tab" aria-selected="true" data-tab="linux">Linux</button>
  <button class="tab" role="tab" aria-selected="false" data-tab="mac">macOS</button>
  <button class="tab" role="tab" aria-selected="false" data-tab="win">Windows</button>
</div>
<div class="tabpanel" data-panel="linux">
<p>Needs a C++20 compiler, CMake ≥ 3.24, Qt ≥ 6.5 (Widgets, Concurrent, Svg), Eigen ≥ 3.4, FFTW 3; GoogleTest for the tests.</p>
${code("shell", `cd software
scripts/build.sh                 # configure + build (Release)
ctest --test-dir build           # 61 tests
build/gui/ssp_gui                # start the interface`, "bash")}</div>
<div class="tabpanel" data-panel="mac" hidden>
<p>Needs the Xcode command-line tools and Homebrew. The script installs the dependencies, builds, bundles Qt into the app, signs it ad hoc and writes a zip. The app has no Apple developer signature: on another Mac use right-click → Open at first launch. <i>This script has not been run on real hardware yet.</i></p>
${code("shell", `xcode-select --install
cd software
bash packaging/macos/build_macos.sh
open packaging/out/ssp-macos-$(uname -m)/ssp_gui.app`, "bash")}</div>
<div class="tabpanel" data-panel="win" hidden>
<p>The tested route is a cross-build in Docker (from Linux, macOS, or Windows with WSL 2). It produces <code>ssp-windows-x64.zip</code> with the programs, all DLLs and the examples — extract and start <code>ssp_gui.exe</code>.</p>
${code("shell", `cd software
packaging/windows/build_windows.sh        # result: packaging/out/ssp-windows-x64.zip`, "bash")}
<p>Native build in an MSYS2 UCRT64 shell (<i>not tested yet</i>):</p>
${code("shell", `pacman -S mingw-w64-ucrt-x86_64-{gcc,cmake,ninja,qt6-base,qt6-svg,eigen3,fftw,pkgconf}
cmake -S software -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DSSP_BUILD_TESTS=OFF
cmake --build build && windeployqt6 build/gui/ssp_gui.exe`, "bash")}</div>
<p class="note">Command line: <code>ssp_aggregate config.json</code> runs a calculation, a scan or a calibration; <code>ssp_aggregate --huckel-tresp molecule.xyz charges.txt 1</code> makes transition charges; <code>ssp_benchmark purcell Au.yml out.csv</code> compares decay rates with exact Mie theory.</p>

<h2 id="validation">Validation</h2>
<p class="lead">Each block is tested against an independent reference. Where the model has a limit, the limit is measured and stated.</p>
<div class="tablewrap"><table><thead><tr><th>Block</th><th>Reference</th><th>Agreement</th></tr></thead><tbody>
<tr><td>Line shapes</td><td>analytic Drude g(t); detailed balance; HEOM for a monomer</td><td>&lt; 1 %</td></tr>
<tr><td>Redfield kinetics</td><td>Boltzmann equilibrium; HEOM for a dimer</td><td>0.03 in population</td></tr>
<tr><td>Modified Redfield</td><td>Redfield limit; detailed balance; HEOM at λ ≈ J</td><td>3 %; 0.5 %; closer than Redfield</td></tr>
<tr><td>Layered sphere</td><td>homogeneous limits to l = 60; Gersten–Nitzan series to l = 90</td><td>10⁻⁶</td></tr>
<tr><td>Aggregate + particle</td><td>transparent particle → bare aggregate spectrum</td><td>&lt; 1 %</td></tr>
<tr><td>Coupling strength g</td><td>analytic Drude sphere</td><td>2 %</td></tr>
<tr><td>Decay near a gold sphere</td><td>exact Mie theory for a dipole</td><td>R ≤ 10 nm: 1–3 % (total), ≤ 10 % (radiative); R = 20 nm: 7–23 %</td></tr>
<tr><td>Fast solver</td><td>exact solver, 300 molecules</td><td>0.2–0.7 %</td></tr>
<tr><td>Monomer fit</td><td>synthetic spectra with known parameters</td><td>within 1–2 σ</td></tr>
<tr><td>C8S3 nanotube J-band</td><td>measured spectrum, no fit to the aggregate</td><td>+100…+200 cm⁻¹ of a −2200 cm⁻¹ shift</td></tr>
</tbody></table></div>
<div class="shots" style="margin-top:14px">
  <figure><img loading="lazy" src="demo/figures/fig6_decay_vs_mie.png" alt="Quantum yield near a gold sphere: model and exact Mie theory"><figcaption>Quantum yield of a molecule near a gold sphere: model (lines) and exact Mie theory (circles).</figcaption></figure>
  <figure><img loading="lazy" src="demo/figures/fig7_heom.png" alt="Populations: HEOM and Redfield"><figcaption>Exciton relaxation in a six-molecule chain: HEOM and Redfield theory.</figcaption></figure>
</div>

<h2 id="roadmap">Roadmap</h2>
<div class="road">
  <div><b>Uncertainty budget</b>Propagate the spread of the inputs (F, n, μ, structure, bath) to positions and widths.</div>
  <div><b>Band widths</b>Periodic boundary conditions along tubes; contributions to the J-band width.</div>
  <div><b>Larger particles</b>Higher-multipole radiation, retardation inside the particle, shells on spheroids.</div>
  <div><b>Vibronic fluorescence</b>Two-particle states in the Holstein model.</div>
  <div><b>Quantum chemistry in the loop</b>A general geometry → TD-DFT → transition-charges pipeline.</div>
  <div><b>More spectroscopies</b>Raman, Brillouin, pump–probe and two-dimensional spectra.</div>
  <div><b>Bayesian calibration</b>Posterior parameters, residual analysis, latent variables.</div>
  <div><b>Atomistic structures</b>Aggregates on particle surfaces from molecular dynamics.</div>
</div>
<footer>SSP Nano guide · generated by <code>docs/guide/build_guide.js</code> · formulas are native MathML (no external scripts). Data used in the examples: Doria et al. 2022, Kuevda et al. 2025, Bondarenko et al. 2020 (all CC BY 4.0); gold optical constants from refractiveindex.info (CC0).</footer>
</main>
</div>
<script>
var MAP = ${JSON.stringify(MAP)};
var LABEL = {done:["Implemented","\\u2713"],partial:["Partly implemented","\\u25D0"],planned:["Not implemented yet","\\u25CB"]};
${fs.readFileSync(path.join(__dirname, "guide_client.js"), "utf8")}
</script>
</body>
</html>
`;
const out = path.join(__dirname, "..", "index.html");
fs.writeFileSync(out, html);
console.log("written", out, (html.length / 1024).toFixed(0) + " kB");
