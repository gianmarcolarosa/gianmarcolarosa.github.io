"use strict";
// Inverse tab: cofactor/adjugate and Gauss-Jordan steps from ../matrix-inverse (math.js lines 1-492, render.js step player).
// Wrapped in a function so its globals (Fr, T, mat, $) do not clash with det.js.
window.INV = (function () {
// Matrix inverse and Moore-Penrose pseudo-inverse: exact rational logic (BigInt fractions).
// Every method returns { steps, ...result }; a step is { text, mats, vars, tone, phase }.

/* ---------- fractions ---------- */
const gcdB = (a, b) => {
  a = a < 0n ? -a : a;
  b = b < 0n ? -b : b;
  while (b) [a, b] = [b, a % b];
  return a;
};
class Fr {
  constructor(n, d = 1n) {
    if (d < 0n) {
      n = -n;
      d = -d;
    }
    const g = gcdB(n, d) || 1n;
    this.n = n / g;
    this.d = d / g;
  }
  add(o) {
    return new Fr(this.n * o.d + o.n * this.d, this.d * o.d);
  }
  sub(o) {
    return new Fr(this.n * o.d - o.n * this.d, this.d * o.d);
  }
  mul(o) {
    return new Fr(this.n * o.n, this.d * o.d);
  }
  div(o) {
    return new Fr(this.n * o.d, this.d * o.n);
  }
  neg() {
    return new Fr(-this.n, this.d);
  }
  abs() {
    return this.n < 0n ? this.neg() : this;
  }
  gt(o) {
    return this.n * o.d > o.n * this.d;
  }
  get zero() {
    return this.n === 0n;
  }
  eq(o) {
    return this.n === o.n && this.d === o.d;
  }
  tex() {
    const sg = this.n < 0n ? "-" : "";
    const a = this.n < 0n ? -this.n : this.n;
    return this.d === 1n ? `${sg}${a}` : `${sg}\\frac{${a}}{${this.d}}`;
  }
  str() {
    return this.d === 1n ? `${this.n}` : `${this.n}/${this.d}`;
  }
  num() {
    // nearest double, safe for large BigInt
    const sh = (x) => {
      const s = x < 0n ? -x : x;
      const bits = s.toString(2).length;
      return bits > 900
        ? [Number(s >> BigInt(bits - 900)), bits - 900]
        : [Number(s), 0];
    };
    const [a, ea] = sh(this.n),
      [b, eb] = sh(this.d);
    return (this.n < 0n ? -1 : 1) * (a / b) * Math.pow(2, ea - eb);
  }
  // accepts "3", "-2.5", ".5", "1/2", "-3/4", "1,5"; otherwise null
  static parse(str) {
    const dec = (t) => {
      const m = /^([+-]?)(\d*)\.?(\d*)$/.exec(t.trim());
      if (!m || (m[2] === "" && m[3] === "")) return null;
      const f = new Fr(
        BigInt((m[2] || "0") + m[3]),
        10n ** BigInt(m[3].length),
      );
      return m[1] === "-" ? f.neg() : f;
    };
    const s = String(str).trim().replace(",", ".");
    if (s.includes("/")) {
      const [a, b, ...rest] = s.split("/");
      if (rest.length) return null;
      const x = dec(a),
        y = dec(b);
      return x && y && !y.zero ? x.div(y) : null;
    }
    return dec(s);
  }
}
const fr = (n, d = 1) => new Fr(BigInt(n), BigInt(d));
const ZERO = fr(0),
  ONE = fr(1);

/* ---------- matrices (arrays of rows of Fr) ---------- */
const mZero = (m, n) =>
  Array.from({ length: m }, () => Array.from({ length: n }, () => ZERO));
const mEye = (n) =>
  mZero(n, n).map((r, i) => r.map((_, j) => (i === j ? ONE : ZERO)));
const mT = (A) => (A.length ? A[0].map((_, j) => A.map((r) => r[j])) : []);
const mMul = (A, B) =>
  A.map((r) =>
    B[0].map((_, j) => r.reduce((s, x, k) => s.add(x.mul(B[k][j])), ZERO)),
  );
const mAdd = (A, B) => A.map((r, i) => r.map((x, j) => x.add(B[i][j])));
const mSub = (A, B) => A.map((r, i) => r.map((x, j) => x.sub(B[i][j])));
const mScale = (A, c) => A.map((r) => r.map((x) => x.mul(c)));
const mEq = (A, B) =>
  A.length === B.length &&
  A.every(
    (r, i) => r.length === B[i].length && r.every((x, j) => x.eq(B[i][j])),
  );
const mIsZero = (A) => A.every((r) => r.every((x) => x.zero));
const mTrace = (A) => A.reduce((s, r, i) => s.add(r[i]), ZERO);
const mCopy = (A) => A.map((r) => r.slice());
const vecCol = (v) => v.map((x) => [x]);
const colOf = (A, j) => A.map((r) => r[j]);
const dot = (u, v) => u.reduce((s, x, i) => s.add(x.mul(v[i])), ZERO);
const minorOf = (A, i, j) =>
  A.filter((_, r) => r !== i).map((r) => r.filter((_, c) => c !== j));
function det(A) {
  const n = A.length;
  if (n === 0) return ONE;
  const W = mCopy(A);
  let d = ONE;
  for (let c = 0; c < n; c++) {
    let p = -1;
    for (let i = c; i < n; i++)
      if (!W[i][c].zero) {
        p = i;
        break;
      }
    if (p < 0) return ZERO;
    if (p !== c) {
      [W[p], W[c]] = [W[c], W[p]];
      d = d.neg();
    }
    d = d.mul(W[c][c]);
    for (let i = c + 1; i < n; i++) {
      if (W[i][c].zero) continue;
      const k = W[i][c].div(W[c][c]);
      for (let j = c; j < n; j++) W[i][j] = W[i][j].sub(k.mul(W[c][j]));
    }
  }
  return d;
}
function rrefQuick(A) {
  const W = mCopy(A),
    m = W.length,
    n = m ? W[0].length : 0,
    piv = [];
  let r = 0;
  for (let c = 0; c < n && r < m; c++) {
    let p = -1;
    for (let i = r; i < m; i++)
      if (!W[i][c].zero) {
        p = i;
        break;
      }
    if (p < 0) continue;
    [W[p], W[r]] = [W[r], W[p]];
    const pv = W[r][c];
    W[r] = W[r].map((x) => x.div(pv));
    for (let i = 0; i < m; i++) {
      if (i === r || W[i][c].zero) continue;
      const k = W[i][c];
      W[i] = W[i].map((x, j) => x.sub(k.mul(W[r][j])));
    }
    piv.push(c);
    r++;
  }
  return { R: W, pivots: piv };
}
const rankOf = (A) => rrefQuick(A).pivots.length;
function invQuick(A) {
  const n = A.length;
  const { R, pivots } = rrefQuick(A.map((r, i) => [...r, ...mEye(n)[i]]));
  if (pivots.length < n || pivots[n - 1] >= n) return null;
  return R.map((r) => r.slice(n));
}

/* ---------- display helpers ---------- */
const fmtNum = (x) => {
  const s = String(+x.toFixed(6));
  return s === "-0" ? "0" : s;
};
function cellTex(x) {
  if (x === null || x === undefined) return "?";
  if (typeof x === "string") return x;
  if (typeof x === "number") return fmtNum(x);
  return x.tex().replace(/\\frac/g, "\\dfrac");
}
// matrix -> LaTeX. hl: {strong:[[i,j]], soft:[[i,j]], dim:[[i,j]], rows:[i], cols:[j]}; split: vertical bar after that column
function matTex(M, o = {}) {
  const hl = o.hl || {};
  const has = (arr, i, j) => (arr || []).some(([a, b]) => a === i && b === j);
  const nc = M.length ? M[0].length : 0;
  const cell = (i, j) => {
    let t = cellTex(M[i][j]);
    // a dimmed cell that is also highlighted keeps the dark text: pale grey on yellow is unreadable
    if (has(hl.dim, i, j) && !has(hl.strong, i, j) && !has(hl.soft, i, j)) t = `\\textcolor{#aab0bd}{${t}}`;
    if (has(hl.strong, i, j)) return `\\colorbox{#ffc24d}{$\\color{#1d2433}{${t}}$}`;
    if (
      has(hl.soft, i, j) ||
      (hl.rows || []).includes(i) ||
      (hl.cols || []).includes(j)
    )
      return `\\colorbox{#fff0b3}{$\\color{#1d2433}{${t}}$}`;
    return t;
  };
  const spec = o.split
    ? "c".repeat(o.split) + "|" + "c".repeat(nc - o.split)
    : "c".repeat(nc);
  const frac = M.some((r) => r.some((x) => /frac/.test(cellTex(x))));
  const body = M.map((r, i) => r.map((_, j) => cell(i, j)).join(" & ")).join(
    frac ? " \\\\[1.3em] " : " \\\\ ",
  );
  return `{\\def\\arraystretch{${frac ? 1.15 : 1.35}}\\left[\\begin{array}{${spec}} ${body} \\end{array}\\right]}`;
}
const T = (s) => `$${s}$`;
const sgn = (f) => (f.n < 0n ? `- ${f.abs().tex()}` : `+ ${f.tex()}`);
const coefTex = (f) =>
  f.eq(ONE) ? "" : f.eq(ONE.neg()) ? "-" : f.tex() + "\\,";
const par = (f) => (f.n < 0n ? `\\left(${f.tex()}\\right)` : f.tex());
const opSwap = (i, j) => `R_{${i + 1}} \\leftrightarrow R_{${j + 1}}`;
const opScale = (i, p) =>
  `R_{${i + 1}} \\leftarrow ${coefTex(ONE.div(p))}R_{${i + 1}}`;
const opElim = (i, j, k) =>
  k.n < 0n
    ? `R_{${i + 1}} \\leftarrow R_{${i + 1}} + ${coefTex(k.abs())}R_{${j + 1}}`
    : `R_{${i + 1}} \\leftarrow R_{${i + 1}} - ${coefTex(k)}R_{${j + 1}}`;
const mat = (name, M, extra = {}) => ({ name, M, ...extra });
const vecInline = (v) =>
  `\\left(${v.map((x) => x.tex()).join(",\\ ")}\\right)^{T}`;
const identityTest = (P) => mEq(P, mEye(P.length));
function sumTex(terms) {
  // terms: strings; joins with + handling a leading minus
  return terms.length ? terms.join(" + ").replace(/\+ -/g, "- ") : "0";
}
function dotTex(u, v) {
  return sumTex(u.map((x, i) => `${par(x)}\\cdot ${par(v[i])}`));
}

/* ---------- final check A * A^-1 = I ---------- */
function verifyStep(A, Inv, nameInv = "A^{-1}") {
  const P = mMul(A, Inv),
    Q = mMul(Inv, A),
    ok = identityTest(P) && identityTest(Q);
  const n = A.length;
  const row = colOf(mT(A), 0);
  const ex = `(A\\,${nameInv})_{11} = ${dotTex(A[0], colOf(Inv, 0))} = ${P[0][0].tex()}`;
  return {
    phase: "Verifica",
    text: `Verifica: moltiplico $A\\,${nameInv}$ riga per colonna con frazioni esatte. Per esempio ${T(ex)}. Il risultato ${ok ? "è la matrice identità $I$: <b>l'inversa è corretta</b> (anche $" + nameInv + "A = I$)" : "<b>NON</b> è l'identità"}.`,
    mats: [
      mat("A", A),
      mat(nameInv, Inv),
      mat("A\\," + nameInv, P, { hl: { strong: P.map((_, i) => [i, i]) } }),
    ],
    vars: { "A\\,A^{-1}": ok ? "I" : "\\ne I" },
    tone: ok ? "ok" : "err",
    ok,
  };
}

/* ---------- Gauss-Jordan engine (also used for rref and for Gram matrices) ---------- */
function gjEngine(W0, o) {
  const W = mCopy(W0),
    m = W.length,
    nc = o.nc;
  const steps = [],
    pivots = [];
  const nm = o.name || "[A\\,|\\,I]";
  const snap = (hl) => mat(nm, mCopy(W), { hl, split: o.split });
  const push = (text, hl, vars, extra = {}) =>
    steps.push({ text, mats: [snap(hl)], vars, phase: o.phase, ...extra });
  let r = 0,
    singular = null;
  for (let col = 0; col < nc && r < m; col++) {
    let p = -1;
    if (o.partial) {
      let best = ZERO;
      for (let i = r; i < m; i++)
        if (W[i][col].abs().gt(best)) {
          best = W[i][col].abs();
          p = i;
        }
    } else {
      for (let i = r; i < m; i++)
        if (!W[i][col].zero) {
          p = i;
          break;
        }
    }
    if (p < 0) {
      if (o.stopOnMissing) {
        singular = { col, r };
        push(
          `Colonna ${col + 1}: dalla riga ${r + 1} in giù tutti gli elementi sono $0$, quindi non c'è nessun <b>pivot</b> (elemento non nullo da usare per azzerare gli altri). Si sono trovati solo ${r} pivot su ${nc}: il rango è minore di ${nc}, la matrice <b>non è invertibile</b> (la colonna ${col + 1} è combinazione delle precedenti).`,
          { cols: [col] },
          { colonna: `${col + 1}`, pivot: "\\text{nessuno}" },
          { tone: "err" },
        );
        break;
      }
      continue;
    }
    if (p !== r) {
      [W[p], W[r]] = [W[r], W[p]];
      push(
        `Colonna ${col + 1}: ${o.partial ? "il pivot scelto è l'elemento di valore assoluto massimo, che sta nella riga " + (p + 1) : "l'elemento in posizione (" + (r + 1) + "," + (col + 1) + ") è $0$: prendo come pivot il primo elemento non nullo, nella riga " + (p + 1)}. Scambio ${T(opSwap(r, p))}; il pivot è ${T(W[r][col].tex())}.`,
        { rows: [r, p], strong: [[r, col]] },
        {
          colonna: `${col + 1}`,
          pivot: W[r][col].tex(),
          operazione: opSwap(r, p),
        },
      );
    }
    const pv = W[r][col];
    if (!pv.eq(ONE)) {
      W[r] = W[r].map((x) => x.div(pv));
      push(
        `Divido la riga ${r + 1} per il pivot ${T(pv.tex())} per ottenere pivot $1$: ${T(opScale(r, pv))}.`,
        { rows: [r], strong: [[r, col]] },
        { colonna: `${col + 1}`, pivot: pv.tex(), operazione: opScale(r, pv) },
      );
    }
    for (let i = 0; i < m; i++) {
      if (i === r || W[i][col].zero) continue;
      const k = W[i][col];
      W[i] = W[i].map((x, j) => x.sub(k.mul(W[r][j])));
      push(
        `Azzero l'elemento ${T(k.tex())} della riga ${i + 1} sotto/sopra il pivot: ${T(opElim(i, r, k))}.`,
        {
          rows: [i],
          cols: [col],
          strong: [
            [r, col],
            [i, col],
          ],
        },
        { colonna: `${col + 1}`, pivot: "1", operazione: opElim(i, r, k) },
      );
    }
    pivots.push(col);
    r++;
  }
  return { steps, W, pivots, singular };
}
const startStep = (name, M, text, split, vars = {}) => ({
  text,
  mats: [mat(name, M, { split })],
  vars: { ...vars },
});

/* ============ PART A, method 1: cofactors / adjugate ============ */
function cofactorSteps(A) {
  const n = A.length,
    steps = [];
  const d = det(A);
  steps.push({
    phase: "Determinante",
    text: `Calcolo il <b>determinante</b> ${T("\\det(A) = " + d.tex())} (con l'eliminazione, esatta). ${d.zero ? "<b>Vale 0: la matrice non è invertibile</b>, perché la formula $A^{-1} = \\operatorname{adj}(A)/\\det(A)$ dividerebbe per zero." : "È diverso da 0: l'inversa esiste e la formula si può applicare."}`,
    mats: [mat("A", A)],
    vars: { "\\det A": d.tex() },
    tone: d.zero ? "err" : undefined,
  });
  if (d.zero) return { steps, inv: null, det: d, reason: "det(A) = 0" };
  const minors = Array.from({ length: n }, () => Array(n).fill(null));
  const cof = Array.from({ length: n }, () => Array(n).fill(null));
  const Mraw = [],
    C = [];
  for (let i = 0; i < n; i++) {
    Mraw.push([]);
    C.push([]);
    for (let j = 0; j < n; j++) {
      const sub = minorOf(A, i, j),
        dm = det(sub),
        s = (i + j) % 2 ? -1 : 1;
      Mraw[i].push(dm);
      C[i].push(s < 0 ? dm.neg() : dm);
    }
  }
  const dimCells = (i, j) => {
    const out = [];
    for (let a = 0; a < n; a++)
      for (let b = 0; b < n; b++) if (a === i || b === j) out.push([a, b]);
    return out;
  };
  const keep = (i, j) => {
    const out = [];
    for (let a = 0; a < n; a++)
      for (let b = 0; b < n; b++) if (a !== i && b !== j) out.push([a, b]);
    return out;
  };
  for (let i = 0; i < n; i++)
    for (let j = 0; j < n; j++) {
      minors[i][j] = Mraw[i][j];
      const subTex = n > 1 ? matTex(minorOf(A, i, j)) : "[\\,]";
      steps.push({
        phase: "Minori",
        text: `Minore ${T(`M_{${i + 1}${j + 1}}`)}: cancello la riga ${i + 1} e la colonna ${j + 1} (in grigio) e calcolo il determinante di ciò che resta: ${T(`M_{${i + 1}${j + 1}} = \\det ${subTex} = ${Mraw[i][j].tex()}`)}. ${n === 1 ? "(Per una matrice $1\\times1$ il minore è il determinante della matrice vuota, cioè $1$.)" : ""} Il cofattore sarà ${T(`C_{${i + 1}${j + 1}} = (-1)^{${i + 1}+${j + 1}} M_{${i + 1}${j + 1}} = ${C[i][j].tex()}`)}.`,
        mats: [
          mat("A", A, {
            hl: { dim: dimCells(i, j), soft: keep(i, j), strong: [[i, j]] },
          }),
          mat("\\text{minori }M", mCopy(minors), { hl: { strong: [[i, j]] } }),
        ],
        vars: {
          "(i,j)": `(${i + 1},${j + 1})`,
          [`M_{ij}`]: Mraw[i][j].tex(),
          [`C_{ij}`]: C[i][j].tex(),
        },
      });
    }
  const S = Array.from({ length: n }, (_, i) =>
    Array.from({ length: n }, (_, j) => ((i + j) % 2 ? "-" : "+")),
  );
  steps.push({
    phase: "Segni",
    text: `I <b>segni</b> $(-1)^{i+j}$ formano una scacchiera: $+$ dove $i+j$ è pari, $-$ dove è dispari.`,
    mats: [mat("(-1)^{i+j}", S), mat("M", Mraw)],
    vars: { "\\text{segni}": "\\text{scacchiera}" },
  });
  steps.push({
    phase: "Cofattori",
    text: `Moltiplico ogni minore per il suo segno: si ottiene la <b>matrice dei cofattori</b> $C$, con ${T("C_{ij} = (-1)^{i+j} M_{ij}")}.`,
    mats: [mat("M", Mraw), mat("C", C)],
    vars: { C: "\\text{cofattori}" },
  });
  const adj = mT(C);
  steps.push({
    phase: "Aggiunta",
    text: `Trasponendo la matrice dei cofattori si ottiene l'<b>aggiunta</b> (o matrice aggiunta classica): ${T("\\operatorname{adj}(A) = C^{T}")}.`,
    mats: [mat("C", C), mat("\\operatorname{adj}(A) = C^{T}", adj)],
    vars: { adj: "C^T" },
  });
  const lap = A[0].map((a, j) => `${par(a)}\\cdot ${par(C[0][j])}`);
  steps.push({
    phase: "Determinante",
    text: `Controllo incrociato: sviluppando il determinante lungo la prima riga (Laplace) ${T("\\det A = \\sum_j a_{1j} C_{1j} = " + sumTex(lap) + " = " + d.tex())}, lo stesso valore di prima.`,
    mats: [
      mat("A", A, { hl: { rows: [0] } }),
      mat("C", C, { hl: { rows: [0] } }),
    ],
    vars: { "\\det A": d.tex() },
  });
  const inv = adj.map((r) => r.map((x) => x.div(d)));
  steps.push({
    phase: "Inversa",
    text: `Divido ogni elemento dell'aggiunta per ${T("\\det A = " + d.tex())}: ${T("A^{-1} = \\dfrac{1}{\\det A}\\operatorname{adj}(A)")}.`,
    mats: [mat("\\operatorname{adj}(A)", adj), mat("A^{-1}", inv)],
    vars: { "A^{-1}": "\\operatorname{adj}(A)/\\det A" },
  });
  steps.push(verifyStep(A, inv));
  return { steps, inv, det: d, cof: C, adj };
}

/* ============ method 2: Gauss-Jordan on [A | I] ============ */
function gaussJordanSteps(A, partial = false, opt = {}) {
  const n = A.length;
  const W0 = A.map((r, i) => [...r, ...mEye(n)[i]]);
  const first = startStep(
    "[A\\,|\\,I]",
    W0,
    `Affianco ad $A$ la matrice identità $I$: scrivo $[A\\,|\\,I]$. Con operazioni sulle righe (scambio due righe, moltiplico una riga per un numero $\\ne0$, sommo a una riga un multiplo di un'altra) porto la parte sinistra a $I$: le stesse operazioni trasformano la parte destra in $A^{-1}$. Pivot scelto: ${partial ? "il più grande in valore assoluto della colonna (<i>pivoting parziale</i>)" : "il primo elemento non nullo"}.`,
    n,
    opt.vars,
  );
  const e = gjEngine(W0, {
    nc: n,
    stopOnMissing: true,
    partial,
    split: n,
    name: "[A\\,|\\,I]",
    phase: opt.phase,
  });
  const steps = [first, ...e.steps];
  if (e.singular) return { steps, inv: null, reason: "rank" };
  const inv = e.W.map((r) => r.slice(n));
  steps.push({
    phase: opt.phase,
    text: `A sinistra c'è $I$: a destra c'è l'inversa $A^{-1}$.`,
    mats: [mat("[I\\,|\\,A^{-1}]", e.W, { split: n }), mat("A^{-1}", inv)],
    vars: { "A^{-1}": "\\text{letta a destra}" },
    tone: "ok",
  });
  if (!opt.noVerify) steps.push(verifyStep(A, inv));
  return { steps, inv };
}

// ---------- rendering helpers: tables, step player, state cards, SVG ----------
const esc = (s) =>
  String(s).replace(
    /[&<>"]/g,
    (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" })[c],
  );
const $ = (id) => document.getElementById(id);
const tex = (s) => "$" + s + "$";
const aligned = (lines) =>
  String.raw`$$\begin{aligned}` +
  lines.join(String.raw` \\ `) +
  String.raw`\end{aligned}$$`;

// Typeset every $...$ and $$...$$ under root with KaTeX (no-op if the CDN is unreachable)
function typeset(root = document.body) {
  if (!window.renderMathInElement) return;
  renderMathInElement(root, {
    strict: "ignore",
    delimiters: [
      { left: "$$", right: "$$", display: true },
      { left: "$", right: "$", display: false },
    ],
  });
  root
    .querySelectorAll(".katex-display")
    .forEach((d) => d.classList.toggle("al", !!d.querySelector(".mtable")));
}

function table(head, rows, cls = "") {
  return (
    `<table class="${cls}"><tr>` +
    head.map((h) => `<th>${h}</th>`).join("") +
    "</tr>" +
    rows
      .map((r) => {
        const cells = r.cells || r;
        return (
          `<tr class="${r.cls || ""}">` +
          cells.map((c) => `<td>${c}</td>`).join("") +
          "</tr>"
        );
      })
      .join("") +
    "</table>"
  );
}

// numbered equations: rows [number, left side, relation + right side, comment]
const numbered = (rows) =>
  String.raw`$$\begin{array}{r l l}` +
  rows
    .map((r) => r.map((c) => String(c).replace(/</g, String.raw`\lt `)))
    .map(
      ([n, l, rr, c]) =>
        String.raw`\textbf{${n}}\ & ${l} & \ ${rr}` +
        (c ? String.raw` \\[-2pt] & {\small ${c}} & ` : ""),
    )
    .join(String.raw` \\[7pt] `) +
  String.raw`\end{array}$$`;

// one matrix box: name above, matrix below (display math)
function matBox(m) {
  const cols = m.M.length ? m.M[0].length : 0;
  const wide = cols > 9 ? "xwide" : cols > 5 ? "wide" : "";
  return `<div class="mbox ${wide}"><div class="mname">${tex(m.name)}</div>$$${matTex(m.M, { hl: m.hl, split: m.split })}$$</div>`;
}

// cards with the variables at step k (changed ones highlighted) + stack of the last states
function stateView(states, cols, k, keep = 8) {
  const cur = states[k],
    prev = k > 0 ? states[k - 1] : null;
  const val = (s, c) => (s[c] === undefined || s[c] === null ? "—" : s[c]);
  const t = (v) => (v === "—" ? "—" : tex(v));
  const lab = (c) => tex(/^[A-Za-z ]{2,}$/.test(c) ? `\\text{${c}}` : c);
  const cards = cols
    .map((c) => {
      const changed = prev && val(prev, c) !== val(cur, c);
      return `<div class="var ${changed ? "changed" : ""}"><div class="vn">${lab(c)}</div><div class="vv">${t(val(cur, c))}</div><div class="vp">${changed && val(prev, c) !== "—" ? "era " + t(val(prev, c)) : "&nbsp;"}</div></div>`;
    })
    .join("");
  const from = Math.max(0, k - keep + 1);
  const frames = states
    .slice(from, k + 1)
    .map((s, i) => ({
      cls: from + i === k ? "cur" : "",
      cells: [
        from + i === 0 ? "inizio" : "passo " + (from + i),
        ...cols.map((c) => t(val(s, c))),
      ],
    }))
    .reverse();
  return `<div class="vars">${cards}</div><h3>Pila degli stati${k + 1 > keep ? ` <span class="muted">(ultimi ${keep})</span>` : ""}</h3><div class="stk">${table(["", ...cols.map((c) => lab(c))], frames)}</div>`;
}

// Step player: slider, prev/next/play, then phase, explanation, matrices, state cards.
// onStep(k, step) lets a section draw something more (e.g. an SVG).
function makeStepper(rootId, delay = 1500, onStep = null) {
  const root = $(rootId);
  root.classList.add("stp");
  const id = (s) => rootId + "_" + s;
  root.innerHTML = `<div class="ctrl"><button id="${id("reset")}" class="alt">⏮</button><button id="${id("prev")}" class="alt">◀</button><button id="${id("play")}">▶ Avvia</button><button id="${id("next")}" class="alt">▶|</button><input id="${id("slider")}" type="range" min="0" max="0" value="0" /><span class="cnt" id="${id("cnt")}"></span></div><div id="${id("phase")}"></div><div class="log" id="${id("log")}"></div><div class="mats" id="${id("mats")}"></div><div id="${id("vars")}"></div>`;
  const c = { steps: [], states: [], cols: [], k: 0, timer: null };
  const el = (s) => $(id(s));
  c.stop = () => {
    clearInterval(c.timer);
    c.timer = null;
    el("play").textContent = "▶ Avvia";
  };
  c.show = (k) => {
    const n = c.steps.length;
    c.k = Math.max(0, Math.min(n - 1, k));
    el("slider").value = c.k;
    const s = c.steps[c.k];
    el("cnt").textContent = `passo ${c.k} di ${n - 1}`;
    el("phase").innerHTML = s.phase
      ? `<span class="phase">${esc(s.phase)}</span>`
      : "";
    el("log").className = "log " + (s.tone || "");
    el("log").innerHTML = s.text;
    el("mats").innerHTML = (s.mats || []).map(matBox).join("");
    const cols = c.cols.filter((x) => c.states[c.k][x] !== undefined);
    el("vars").innerHTML = cols.length ? stateView(c.states, cols, c.k) : "";
    typeset(root);
    if (onStep) onStep(c.k, s, c);
    if (c.k >= n - 1) c.stop();
  };
  c.set = (steps) => {
    c.stop();
    c.steps = steps;
    const cols = [];
    c.states = [];
    let st = {};
    steps.forEach((s) => {
      st = { ...st, ...(s.vars || {}) };
      c.states.push(st);
      Object.keys(s.vars || {}).forEach(
        (k) => cols.includes(k) || cols.push(k),
      );
    });
    c.cols = cols;
    el("slider").max = steps.length - 1;
    c.show(0);
  };
  el("slider").addEventListener("input", () => {
    c.stop();
    c.show(+el("slider").value);
  });
  el("prev").addEventListener("click", () => {
    c.stop();
    c.show(c.k - 1);
  });
  el("next").addEventListener("click", () => {
    c.stop();
    c.show(c.k + 1);
  });
  el("reset").addEventListener("click", () => {
    c.stop();
    c.show(0);
  });
  el("play").addEventListener("click", () => {
    if (c.timer) return c.stop();
    if (c.k >= c.steps.length - 1) c.show(0);
    el("play").textContent = "⏸ Pausa";
    c.timer = setInterval(() => c.show(c.k + 1), delay);
  });
  return c;
}

return { Fr, fr, matTex, cofactorSteps, gaussJordanSteps, makeStepper, typeset, numbered };
})();
