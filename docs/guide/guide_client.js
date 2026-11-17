// Client script of the guide: program map, filters, tabs, theme switch and three small calculators.
(function () {
  var $ = function (s, r) { return (r || document).querySelector(s); };
  var el = function (tag, cls, text) { var e = document.createElement(tag); if (cls) e.className = cls; if (text != null) e.textContent = text; return e; };
  var css = function (v) { return getComputedStyle(document.documentElement).getPropertyValue(v).trim(); };

  // ---------------- program map
  var grid = $("#mapgrid"), detail = $("#detail"), counts = { all: 0, done: 0, partial: 0, planned: 0 }, current = null;
  MAP.forEach(function (st) {
    var box = el("section", "stage"), h = el("h4", null, st.stage), n = el("span");
    h.appendChild(n); box.appendChild(h);
    var doneN = 0;
    st.items.forEach(function (it) {
      var b = el("button", "node " + it[1]);
      b.type = "button"; b.dataset.status = it[1]; b.setAttribute("aria-expanded", "false");
      var dot = el("span", "dot", LABEL[it[1]][1]); dot.title = LABEL[it[1]][0];
      b.appendChild(dot); b.appendChild(el("span", null, it[0]));
      b.addEventListener("click", function () { show(it, b); });
      box.appendChild(b);
      counts.all++; counts[it[1]]++; if (it[1] === "done") doneN++;
    });
    n.textContent = doneN + "/" + st.items.length;
    n.title = "implemented / total";
    grid.appendChild(box);
  });
  ["all", "done", "partial", "planned"].forEach(function (k) { $("#c-" + k).textContent = "(" + counts[k] + ")"; });

  function show(it, btn) {
    if (current) current.setAttribute("aria-expanded", "false");
    current = btn; btn.setAttribute("aria-expanded", "true");
    detail.textContent = "";
    var h = el("h4", null, it[0]), bd = el("span", "badge " + it[1]);
    bd.appendChild(el("i", null, LABEL[it[1]][1])); bd.appendChild(document.createTextNode(LABEL[it[1]][0]));
    h.appendChild(bd); detail.appendChild(h);
    var dl = el("dl");
    dl.appendChild(el("dt", null, "Method")); dl.appendChild(el("dd", null, it[2]));
    dl.appendChild(el("dt", null, "Gives")); dl.appendChild(el("dd", null, it[3]));
    if (it[4]) {
      dl.appendChild(el("dt", null, "Details"));
      var dd = el("dd"), a = el("a", null, "formulas and code →");
      a.href = (it[4] === "config" || it[4] === "geometry") ? "#config" : "#m-" + it[4];
      dd.appendChild(a); dl.appendChild(dd);
    }
    detail.appendChild(dl);
  }
  document.querySelectorAll(".chip").forEach(function (c) {
    c.addEventListener("click", function () {
      document.querySelectorAll(".chip").forEach(function (x) { x.setAttribute("aria-pressed", String(x === c)); });
      var f = c.dataset.f;
      document.querySelectorAll(".node").forEach(function (n) { n.hidden = !(f === "all" || n.dataset.status === f); });
    });
  });

  // ---------------- tabs
  document.querySelectorAll(".tab").forEach(function (t) {
    t.addEventListener("click", function () {
      document.querySelectorAll(".tab").forEach(function (x) { x.setAttribute("aria-selected", String(x === t)); });
      document.querySelectorAll(".tabpanel").forEach(function (p) { p.hidden = p.dataset.panel !== t.dataset.tab; });
    });
  });

  // ---------------- theme
  var redraws = [];
  $("#theme").addEventListener("click", function () {
    var root = document.documentElement, dark = root.dataset.theme ? root.dataset.theme === "dark" : matchMedia("(prefers-color-scheme: dark)").matches;
    root.dataset.theme = dark ? "light" : "dark";
    redraws.forEach(function (f) { f(); });
  });

  // ---------------- calculators
  function widget(slot, title, sliders, draw) {
    var w = el("div", "widget"), left = el("div"), right = el("div");
    w.appendChild(el("h5", null, title));
    var vals = {}, outs = {};
    sliders.forEach(function (s) {
      var lab = el("label"), out = el("output"), inp = el("input");
      lab.appendChild(el("span", null, s.label)); lab.appendChild(out);
      inp.type = "range"; inp.min = s.min; inp.max = s.max; inp.step = s.step; inp.value = s.value;
      inp.setAttribute("aria-label", s.label);
      lab.appendChild(inp); left.appendChild(lab);
      vals[s.key] = s.value; outs[s.key] = [out, s];
      inp.addEventListener("input", function () { vals[s.key] = parseFloat(inp.value); update(); });
    });
    var read = el("div", "readout"); left.appendChild(read);
    var cv = el("canvas"); cv.setAttribute("role", "img"); cv.setAttribute("aria-label", title); right.appendChild(cv);
    w.appendChild(left); w.appendChild(right); slot.appendChild(w);
    function update() {
      Object.keys(outs).forEach(function (k) { outs[k][0].textContent = vals[k].toFixed(outs[k][1].digits) + " " + outs[k][1].unit; });
      var r = cv.getBoundingClientRect(), d = window.devicePixelRatio || 1;
      cv.width = Math.max(200, r.width * d); cv.height = Math.max(120, r.height * d);
      var g = cv.getContext("2d"); g.setTransform(d, 0, 0, d, 0, 0); g.clearRect(0, 0, r.width, r.height);
      g.font = "12px system-ui, sans-serif";
      draw(vals, g, r.width, r.height, read, { ink: css("--ink"), ink2: css("--ink2"), line: css("--line"), a: css("--accent"), b: css("--accent2"), ok: css("--done") });
    }
    redraws.push(update); window.addEventListener("resize", update); update();
  }
  function axes(g, W, H, C, x0, x1, y0, y1, xl, yl) {
    var L = 66, R = 16, T = 12, B = 34;
    var X = function (x) { return L + (x - x0) / (x1 - x0) * (W - L - R); }, Y = function (y) { return H - B - (y - y0) / (y1 - y0) * (H - T - B); };
    g.strokeStyle = C.line; g.lineWidth = 1; g.strokeRect(L, T, W - L - R, H - T - B);
    g.fillStyle = C.ink2; g.textAlign = "center"; g.fillText(xl, (L + W - R) / 2, H - 6);
    g.save(); g.translate(12, (T + H - B) / 2); g.rotate(-Math.PI / 2); g.fillText(yl, 0, 0); g.restore();
    for (var i = 0; i <= 4; i++) {
      var xv = x0 + (x1 - x0) * i / 4, yv = y0 + (y1 - y0) * i / 4;
      g.textAlign = i === 0 ? "left" : i === 4 ? "right" : "center"; g.fillText(Math.round(xv), X(xv), H - B + 14);
      g.textAlign = "right"; g.fillText(Math.round(yv), L - 5, Y(yv) + 4);
    }
    return { X: X, Y: Y };
  }
  function curve(g, P, f, x0, x1, col, dash) {
    g.beginPath(); g.setLineDash(dash || []);
    for (var i = 0; i <= 200; i++) { var x = x0 + (x1 - x0) * i / 200, y = f(x); if (i) g.lineTo(P.X(x), P.Y(y)); else g.moveTo(P.X(x), P.Y(y)); }
    g.strokeStyle = col; g.lineWidth = 2; g.stroke(); g.setLineDash([]);
  }

  // 1. point-dipole coupling of two parallel dipoles
  var s1 = $('[data-widget="w-coupling"]');
  if (s1) widget(s1, "Try it: coupling of two parallel transition dipoles", [
    { key: "r", label: "distance r", min: 0.4, max: 3, step: 0.05, value: 1.0, digits: 2, unit: "nm" },
    { key: "th", label: "angle θ between dipoles and the line joining them", min: 0, max: 90, step: 1, value: 20, digits: 0, unit: "°" },
    { key: "mu", label: "transition dipole μ", min: 1, max: 20, step: 0.5, value: 10, digits: 1, unit: "D" },
    { key: "fs", label: "screening f_s", min: 0.2, max: 1, step: 0.02, value: 1, digits: 2, unit: "" }],
  function (v, g, W, H, read, C) {
    var J = function (th) { var c = Math.cos(th * Math.PI / 180); return v.fs * 5.034 * v.mu * v.mu * (1 - 3 * c * c) / Math.pow(v.r, 3); };
    var m = Math.max(Math.abs(J(0)), Math.abs(J(90)), 1);
    var P = axes(g, W, H, C, 0, 90, -m, m, "θ, degrees", "J, cm⁻¹");
    g.strokeStyle = C.ink2; g.setLineDash([3, 3]); g.beginPath(); g.moveTo(P.X(0), P.Y(0)); g.lineTo(P.X(90), P.Y(0)); g.moveTo(P.X(54.7), P.Y(-m)); g.lineTo(P.X(54.7), P.Y(m)); g.stroke(); g.setLineDash([]);
    g.fillStyle = C.ink2; g.textAlign = "left"; g.fillText("magic angle 54.7°", P.X(54.7) + 4, P.Y(m) + 14);
    curve(g, P, J, 0, 90, C.a);
    var j = J(v.th); g.fillStyle = C.b; g.beginPath(); g.arc(P.X(v.th), P.Y(j), 5, 0, 7); g.fill();
    read.innerHTML = "J = <b>" + j.toFixed(0) + " cm⁻¹</b> — " + (j < 0 ? "J-type: head-to-tail, red-shifted bright band" : "H-type: side-by-side, blue-shifted bright band");
  });

  // 2. coupled oscillators
  var s2 = $('[data-widget="w-anticross"]');
  if (s2) widget(s2, "Try it: plasmon–exciton anticrossing and the strong-coupling criteria", [
    { key: "g", label: "coupling g", min: 0, max: 1500, step: 10, value: 300, digits: 0, unit: "cm⁻¹" },
    { key: "gp", label: "plasmon width γp (FWHM)", min: 200, max: 3000, step: 20, value: 1600, digits: 0, unit: "cm⁻¹" },
    { key: "gx", label: "exciton width γx (FWHM)", min: 50, max: 1500, step: 10, value: 400, digits: 0, unit: "cm⁻¹" }],
  function (v, g, W, H, read, C) {
    var wp = 19200, x0 = 17000, x1 = 21400;
    var br = function (wx, sgn) {
      var hr = 0.5 * (wp - wx), hi = -0.25 * (v.gp - v.gx);
      var re = v.g * v.g + hr * hr - hi * hi, imv = 2 * hr * hi, mod = Math.sqrt(re * re + imv * imv);
      var sr = Math.sqrt(Math.max(0, (mod + re) / 2));
      return 0.5 * (wp + wx) + sgn * sr;
    };
    var P = axes(g, W, H, C, x0, x1, 16400, 22000, "exciton energy, cm⁻¹", "hybrid modes, cm⁻¹");
    curve(g, P, function () { return wp; }, x0, x1, C.ink2, [4, 4]);
    curve(g, P, function (x) { return x; }, x0, x1, C.ink2, [4, 4]);
    curve(g, P, function (x) { return br(x, 1); }, x0, x1, C.a);
    curve(g, P, function (x) { return br(x, -1); }, x0, x1, C.b);
    var vis = 2 * v.g > 0.5 * (v.gp + v.gx), anti = v.g > Math.abs(v.gp - v.gx) / 4, coop = 4 * v.g * v.g / (v.gp * v.gx);
    read.innerHTML = "2g = <b>" + (2 * v.g).toFixed(0) + "</b> vs (γp+γx)/2 = <b>" + (0.5 * (v.gp + v.gx)).toFixed(0) + "</b> → splitting " + (vis ? "<b>visible</b>" : "not visible") +
      "<br>anticrossing: <b>" + (anti ? "yes" : "no") + "</b> · cooperativity C = <b>" + coop.toFixed(2) + "</b> → " + (vis ? "strong coupling" : anti ? "intermediate coupling" : "weak coupling");
  });

  // 3. dispersive site energy
  var s3 = $('[data-widget="w-shift"]');
  if (s3) widget(s3, "Try it: where does a molecule absorb inside an aggregate?", [
    { key: "F", label: "Bayliss factor F", min: 0.5, max: 2.5, step: 0.01, value: 1.43, digits: 2, unit: "eV" },
    { key: "na", label: "refractive index of the dye layer", min: 1.3, max: 2.2, step: 0.01, value: 1.78, digits: 2, unit: "" },
    { key: "ns", label: "refractive index of the monomer's solvent", min: 1.0, max: 1.7, step: 0.01, value: 1.33, digits: 2, unit: "" },
    { key: "s", label: "local packing ΣV² / ⟨ΣV²⟩ (1 = continuum)", min: 0.3, max: 1.6, step: 0.01, value: 1.0, digits: 2, unit: "" }],
  function (v, g, W, H, read, C) {
    var f = function (n) { return (n * n - 1) / (2 * n * n + 1); }, E0 = 19900, eV = 8065.544;
    var gas = E0 + v.F * f(v.ns) * eV, agg = gas - v.s * v.F * f(v.na) * eV;
    var lo = 15000, hi = 23500, L = 60, Rr = 20, Y = function (e) { return H - 28 - (e - lo) / (hi - lo) * (H - 50); };
    var cols = [["monomer in solvent", E0, C.ink2], ["gas phase", gas, C.a], ["in the aggregate", agg, C.b]], cw = (W - L - Rr) / 3;
    g.strokeStyle = C.line; g.beginPath(); g.moveTo(L - 10, Y(lo)); g.lineTo(L - 10, Y(hi)); g.stroke();
    g.fillStyle = C.ink2; g.textAlign = "right";
    for (var e = 16000; e <= 23000; e += 1000) g.fillText(e, L - 14, Y(e) + 4);
    cols.forEach(function (c, i) {
      var x = L + i * cw + 10, y = Y(Math.min(hi, Math.max(lo, c[1])));
      g.strokeStyle = c[2]; g.lineWidth = 3; g.beginPath(); g.moveTo(x, y); g.lineTo(x + cw - 20, y); g.stroke();
      g.fillStyle = C.ink; g.textAlign = "center"; g.fillText(Math.round(c[1]) + " cm⁻¹", x + (cw - 20) / 2, y - 7);
      g.fillStyle = C.ink2; g.fillText(c[0], x + (cw - 20) / 2, H - 8);
    });
    read.innerHTML = "back to the gas phase: <b>+" + (gas - E0).toFixed(0) + " cm⁻¹</b> · dispersive shift: <b>" + (agg - gas).toFixed(0) + " cm⁻¹</b><br>net shift of the site energy: <b>" + (agg - E0).toFixed(0) + " cm⁻¹</b> → " + (1e7 / agg).toFixed(0) + " nm";
  });
})();
