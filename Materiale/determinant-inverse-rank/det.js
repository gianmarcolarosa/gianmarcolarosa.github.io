"use strict";
// Determinant logic: exact fraction arithmetic (BigInt)
// and HTML generation for every step. No DOM dependency.

/* ---------- Fractions ---------- */
const gcd = (a, b) => {
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
    const g = gcd(n, d) || 1n;
    this.n = n / g;
    this.d = d / g;
  }
  add(o) {
    return new Fr(this.n * o.d + o.n * this.d, this.d * o.d);
  }
  sub(o) {
    return this.add(o.neg());
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
  html() {
    return m(this.tex());
  }
  dec() {
    const v = Number(this.n) / Number(this.d);
    return String(+v.toFixed(6));
  }
  // accepts "3", "-2.5", ".5", "1/2", "-3/4"; otherwise null
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
const ONE = new Fr(1n),
  ZERO = new Fr(0n);

/* ---------- HTML helpers ---------- */
const clone = (M) => M.map((r) => r.slice());
const m = (t) => `\\(${t}\\)`; // inline LaTeX (rendered by KaTeX)
const T = (x) => x.tex();
const PT = (x) => (x.n < 0n ? `(${x.tex()})` : x.tex()); // parentheses if negative
const FAC_COLORS = [
  "#e8590c",
  "#1c7ed6",
  "#2f9e44",
  "#ae3ec9",
  "#d6336c",
  "#0c8599",
];
const R = (i) => `R_{${i + 1}}`;

function mat(M, cls = () => "") {
  const rows = M.map(
    (row, r) =>
      `<tr>${row.map((x, c) => `<td class="${cls(r, c)}">${x.html()}</td>`).join("")}</tr>`,
  );
  const size = M.length >= 6 ? "xs" : M.length >= 4 ? "sm" : ""; // smaller cells for big matrices
  return `<div class="matw"><table class="mat ${size}">${rows.join("")}</table></div>`;
}

function step(title, formula, body, note = "", cls = "") {
  return (
    `<div class="step ${cls}"><div class="stitle">${title}</div>` +
    (formula ? `<div class="sform">${formula}</div>` : "") +
    `<div class="sbody">${body}</div>` +
    (note ? `<div class="snote">${note}</div>` : "") +
    `</div>`
  );
}

/* ---------- Gauss and Gauss-Jordan ---------- */
function gauss(A, jordan, partial = false) {
  const n = A.length,
    M = clone(A),
    out = [];
  out.push(
    step(
      "Matrice iniziale",
      jordan
        ? "Riduco A alla matrice identità, tenendo traccia dei fattori che moltiplicano il determinante."
        : "Riduco A a forma triangolare superiore con operazioni elementari di riga.",
      mat(M),
    ),
  );
  let sw = 0,
    acc = ONE;
  const facs = []; // determinant factors (Gauss-Jordan)
  let count = 1;

  for (let k = 0; k < n; k++) {
    let r = k;
    if (partial) {
      for (let i = k + 1; i < n; i++)
        if (M[i][k].abs().gt(M[r][k].abs())) r = i;
    } else while (r < n && M[r][k].zero) r++;
    if (r === n || M[r][k].zero) {
      out.push(
        step(
          `Passo ${count++}: colonna ${k + 1} nulla`,
          `Dalla riga ${k + 1} in giù la colonna ${k + 1} contiene solo zeri.`,
          mat(M, (i, j) => (j === k && i >= k ? "chg" : "")),
          `Non esiste pivot ⇒ le righe sono linearmente dipendenti ⇒ <b>${m("\\det(A) = 0")}</b>.`,
        ),
      );
      return { det: ZERO, html: out.join("") };
    }
    if (r !== k) {
      const why = M[k][k].zero
        ? `Pivot nullo in posizione (${k + 1},${k + 1}): scambio le righe.`
        : `Pivoting parziale: nella colonna ${k + 1} il pivot di modulo maggiore (${m(PT(M[r][k]))}) è nella riga ${r + 1}, lo porto in posizione (${k + 1},${k + 1}).`;
      [M[k], M[r]] = [M[r], M[k]];
      sw++;
      if (jordan) {
        acc = acc.neg();
        facs.push({
          f: "(-1)",
          why: `scambio ${m(`${R(k)} \\leftrightarrow ${R(r)}`)} (cambia segno)`,
          n: count,
        });
      }
      out.push(
        step(
          `Passo ${count++}: scambio di righe`,
          m(`${R(k)} \\leftrightarrow ${R(r)}`),
          mat(M, (i) => (i === k || i === r ? "chg" : "")),
          `${why} <b>Lo scambio cambia il segno del determinante</b> (scambi finora: ${sw}).`,
        ),
      );
    }
    if (jordan) {
      const p = M[k][k];
      if (!p.eq(ONE)) {
        M[k] = M[k].map((x) => x.div(p));
        acc = acc.mul(p);
        facs.push({
          f: PT(p),
          why: `${m(`${R(k)} \\leftarrow \\frac{${R(k)}}{${T(p)}}`)} (divisione per il pivot)`,
          n: count,
        });
        out.push(
          step(
            `Passo ${count++}: normalizzo il pivot`,
            m(`${R(k)} \\leftarrow \\frac{${R(k)}}{${T(p)}}`),
            mat(M, (i, j) => (i === k ? (j === k ? "pivot" : "chg") : "")),
            `Dividere una riga per ${m(PT(p))} moltiplica il determinante per ${m(`\\frac{1}{${T(p)}}`)}, quindi tengo il fattore: ${m(`\\det(A) = ${PT(p)} \\cdot \\det(A')`)}. Fattore accumulato: ${m(T(acc))}.`,
          ),
        );
      }
    }
    for (let i = jordan ? 0 : k + 1; i < n; i++) {
      if (i === k || M[i][k].zero) continue;
      const mult = M[i][k].div(M[k][k]);
      M[i] = M[i].map((x, j) => x.sub(mult.mul(M[k][j])));
      out.push(
        step(
          `Passo ${count++}: azzero la posizione (${i + 1},${k + 1})`,
          m(`${R(i)} \\leftarrow ${R(i)} - ${PT(mult)} \\cdot ${R(k)}`),
          mat(M, (a, b) =>
            a === i ? (b === k ? "zero" : "chg") : a === k ? "pivot" : "",
          ),
          "<b>Questa operazione non cambia il determinante.</b>",
        ),
      );
    }
  }

  if (jordan) {
    const col = (i) => FAC_COLORS[i % FAC_COLORS.length];
    const paint = (x, i) => `\\textcolor{${col(i)}}{${x.f}}`;
    const f = facs.length ? facs.map(paint).join(" \\cdot ") : "1";
    const cards = facs.length
      ? facs
          .map(
            (x, i) =>
              `<div class="fac" style="border-color:${col(i)}"><div class="fn">Passo ${x.n}</div>` +
              `<div>${x.why}</div><div class="fv">fattore ${m(paint(x, i))}</div></div>`,
          )
          .join("")
      : '<div class="fac"><div>Nessun fattore: la matrice era già l\'identità.</div></div>';
    out.push(
      step(
        "Matrice identità raggiunta",
        m(`\\displaystyle \\det(A) = ${f} = ${T(acc)}`),
        `<div class="jfinal">${mat(M, (i, j) => (i === j ? "pivot" : ""))}<div class="facs">${cards}</div></div>`,
        `${m("\\det(I) = 1")}, quindi ${m("\\det(A)")} è il prodotto dei <b>fattori colorati</b>: uno per ogni scambio di righe ${m("(-1)")} e uno per ogni divisione di riga per il pivot. Le sottrazioni di righe non cambiano il determinante e non contano.`,
        "wide",
      ),
    );
    return { det: acc, html: out.join("") };
  }
  const diag = M.map((r, i) => r[i]);
  const prod = diag.reduce((a, b) => a.mul(b), ONE);
  const det = sw % 2 ? prod.neg() : prod;
  const B = "#1c7ed6",
    O = "#e8590c";
  const col = (c, t) => `\\textcolor{${c}}{${t}}`;
  const sgn = sw % 2 ? "-1" : "1";
  const card = (c, title, tex, extra = "") =>
    `<div class="fac" style="border-color:${c}"><div class="fn">${title}</div>` +
    `<div class="fv">${m(`\\displaystyle ${col(c, tex)}`)}</div>${extra ? `<div class="fn">${extra}</div>` : ""}</div>`;
  out.push(
    step(
      "Matrice triangolare superiore",
      m(
        `\\displaystyle \\det(A) = ${col(O, `(-1)^{${sw}}`)} \\cdot ${col(B, diag.map(PT).join(" \\cdot "))} = ${col(O, sgn)} \\cdot ${col(B, PT(prod))} = ${T(det)}`,
      ),
      `<div class="jfinal">${mat(M, (i, j) => (i === j ? "pivot" : ""))}<div class="facs">` +
        card(
          B,
          "Prodotto degli elementi della diagonale",
          `${diag.map(PT).join(" \\cdot ")} = ${T(prod)}`,
          "In una matrice triangolare il determinante è il prodotto della diagonale (celle azzurre).",
        ) +
        card(
          O,
          `Segno: ${sw} scambi di righe`,
          `(-1)^{${sw}} = ${sgn}`,
          sw === 0
            ? "Nessuno scambio: il segno non cambia, perché l'esponente è il numero di scambi."
            : `Ogni scambio cambia il segno: ${sw % 2 ? "numero dispari, il segno si inverte" : "numero pari, i cambi si annullano a coppie"}.`,
        ) +
        `</div></div>`,
      `Le sottrazioni di righe non cambiano il determinante: contano solo gli scambi (segno) e la diagonale finale.`,
      "wide",
    ),
  );
  return { det, html: out.join("") };
}

/* ---------- Sarrus (3×3) ---------- */
function sarrus(M) {
  const out = [];
  const ext = M.map((r) => [...r, r[0], r[1]]);
  const ghost = (c) => (c >= 3 ? "ghost " : "");
  out.push(
    step(
      "Affianco le prime due colonne",
      "Riscrivo a destra della matrice la 1ª e la 2ª colonna.",
      mat(ext, (r, c) => ghost(c).trim()),
    ),
  );

  const pos = [0, 1, 2].map((i) => [0, 1, 2].map((r) => ext[r][r + i]));
  const neg = [0, 1, 2].map((i) => [0, 1, 2].map((r) => ext[r][i + 2 - r]));
  const prod = (v) => v.reduce((a, b) => a.mul(b), ONE);
  const pp = pos.map(prod),
    np = neg.map(prod);

  out.push(
    step(
      "Diagonali discendenti (segno +)",
      pos
        .map(
          (v, i) =>
            `<span class="t p${i}">${m(`${v.map(PT).join(" \\cdot ")} = ${T(pp[i])}`)}</span>`,
        )
        .join(" &nbsp; "),
      mat(ext, (r, c) => {
        for (let i = 0; i < 3; i++)
          if (c === r + i) return `${ghost(c)}p${i}`.trim();
        return ghost(c).trim();
      }),
      "Ogni diagonale che scende da sinistra a destra dà un prodotto con segno positivo.",
    ),
  );

  out.push(
    step(
      "Diagonali ascendenti (segno −)",
      neg
        .map(
          (v, i) =>
            `<span class="t n${i}">${m(`${v.map(PT).join(" \\cdot ")} = ${T(np[i])}`)}</span>`,
        )
        .join(" &nbsp; "),
      mat(ext, (r, c) => {
        for (let i = 0; i < 3; i++)
          if (c === i + 2 - r) return `${ghost(c)}n${i}`.trim();
        return ghost(c).trim();
      }),
      "Ogni diagonale che sale da sinistra a destra dà un prodotto con segno negativo.",
    ),
  );

  const sp = pp.reduce((a, b) => a.add(b), ZERO),
    sn = np.reduce((a, b) => a.add(b), ZERO);
  const det = sp.sub(sn);
  const G = "#2f9e44",
    Rd = "#e03131";
  const col = (c, t) => `\\textcolor{${c}}{${t}}`;
  const card = (c, title, terms, sum) =>
    `<div class="fac" style="border-color:${c}"><div class="fn">${title}</div>` +
    `<div class="fv">${m(`\\displaystyle ${col(c, `${terms.map(PT).join(" + ")} = ${T(sum)}`)}`)}</div></div>`;
  out.push(
    step(
      "Somma finale",
      m(
        `\\displaystyle \\det(A) = ${col(G, `\\left(${pp.map(PT).join(" + ")}\\right)`)} - ${col(Rd, `\\left(${np.map(PT).join(" + ")}\\right)`)}`,
      ) +
        "<br>" +
        m(
          `\\displaystyle \\det(A) = ${col(G, T(sp))} - ${col(Rd, PT(sn))} = ${T(det)}`,
        ),
      `<div class="jfinal">${mat(M)}<div class="facs">` +
        card(G, "Somma delle diagonali discendenti (+)", pp, sp) +
        card(Rd, "Somma delle diagonali ascendenti (−)", np, sn) +
        `</div></div>`,
      `Il determinante è la somma dei tre prodotti <b style="color:${G}">verdi</b> meno la somma dei tre prodotti <b style="color:${Rd}">rossi</b>.`,
      "wide",
    ),
  );
  return { det, html: out.join("") };
}

/* ---------- 2×2 formula ---------- */
function det2(M) {
  const [[a, b], [c, d]] = M;
  const ad = a.mul(d),
    bc = b.mul(c),
    det = ad.sub(bc);
  const html =
    step(
      "Formula 2×2",
      m("\\det\\begin{pmatrix} a & b \\\\ c & d \\end{pmatrix} = a d - b c"),
      mat(M, (r, k) => (r === k ? "p0" : "n0")),
      "Diagonale principale (verde) meno diagonale secondaria (rosso).",
    ) +
    step(
      "Sostituisco i valori",
      m(
        `\\det(A) = ${PT(a)} \\cdot ${PT(d)} - ${PT(b)} \\cdot ${PT(c)} = ${T(ad)} - ${PT(bc)} = ${T(det)}`,
      ),
      "",
    );
  return { det, html };
}

/* ---------- Laplace ---------- */
const minor = (M, i, j) =>
  M.filter((_, r) => r !== i).map((row) => row.filter((_, c) => c !== j));

function bestLine(M) {
  const n = M.length;
  let best = { t: "row", i: 0, z: -1 };
  for (let i = 0; i < n; i++) {
    const z = M[i].filter((x) => x.zero).length;
    if (z > best.z) best = { t: "row", i, z };
  }
  for (let j = 0; j < n; j++) {
    const z = M.filter((r) => r[j].zero).length;
    if (z > best.z) best = { t: "col", i: j, z };
  }
  return best;
}

// sign pattern (-1)^(i+j), the chosen line highlighted
function signs(n, isRow, L) {
  const rows = Array.from(
    { length: n },
    (_, i) =>
      `<tr>${Array.from(
        { length: n },
        (_, j) =>
          `<td class="${(isRow ? i === L : j === L) ? "hl" : ""}">${(i + j) % 2 ? "−" : "+"}</td>`,
      ).join("")}</tr>`,
  );
  return `<div class="matw signs"><table class="mat">${rows.join("")}</table></div>`;
}

// pick = {t:'row'|'col', i}; top = first level (does not shortcut 2×2); depth = nesting level
function laplace(M, pick, top = true, depth = 0) {
  const n = M.length;
  if (n === 1) return { det: M[0][0], html: "" };
  if (n === 2 && !top) {
    const [[a, b], [c, d]] = M;
    const det = a.mul(d).sub(b.mul(c));
    return {
      det,
      html: `<div class="tres">${m(`\\det = ${PT(a)} \\cdot ${PT(d)} - ${PT(b)} \\cdot ${PT(c)} = ${T(det)}`)}</div>`,
    };
  }
  const auto = !pick;
  pick = pick || bestLine(M);
  const isRow = pick.t === "row",
    L = pick.i;
  const name = isRow ? "riga" : "colonna";
  const cellsOf = (k) => (isRow ? [L, k] : [k, L]);
  const parts = [];
  parts.push(
    `<div class="lhead">` +
      (depth
        ? `Sviluppo lungo la <b>${name} ${L + 1}</b> (${auto ? "la più ricca di zeri" : "scelta"}).`
        : `<div class="lformula">${m(`\\displaystyle \\det(A) = \\sum_{${isRow ? "j" : "i"}} (-1)^{i+j} \\, a_{ij} \\, \\det(M_{ij})`)}</div>` +
          `Sviluppo lungo la <b>${name} ${L + 1}</b>${auto ? " (la più ricca di zeri)" : ""}. ` +
          `${m("M_{ij}")} è il <b>minore</b>: la matrice senza la riga ${m("i")} e la colonna ${m("j")}. ` +
          `Il segno ${m("(-1)^{i+j}")} segue la scacchiera qui a fianco. I termini con ${m("a_{ij} = 0")} sono nulli e si omettono.`) +
      `</div>`,
  );
  parts.push(
    `<div class="trow">${mat(M, (r, c) => ((isRow ? r === L : c === L) ? "hl" : ""))}` +
      `<div class="slabel">segni<br>${signs(n, isRow, L)}</div></div>`,
  );

  const terms = [];
  let total = ZERO;
  let idx = 0;
  for (let k = 0; k < n; k++) {
    const [i, j] = cellsOf(k);
    const a = M[i][j];
    const s = (i + j) % 2 === 0 ? "+" : "-";
    const id = `${i + 1}${j + 1}`;
    if (a.zero) {
      parts.push(
        `<div class="skip">${m(`a_{${id}} = 0`)} ⇒ termine nullo, omesso.</div>`,
      );
      continue;
    }
    const color = FAC_COLORS[idx++ % FAC_COLORS.length];
    const sub_ = minor(M, i, j);
    const r = laplace(sub_, null, false, depth + 1);
    const val = r.det;
    const contrib = s === "+" ? a.mul(val) : a.mul(val).neg();
    terms.push({ s, a, val, color });
    total = total.add(contrib);
    const sz = sub_.length;
    parts.push(
      `<details class="term" style="--tc:${color}"${depth === 0 ? " open" : ""}>` +
        `<summary><span class="thead">${m(`\\textcolor{${color}}{${s}\\, a_{${id}} \\cdot \\det(M_{${id}})}`)}` +
        `<span class="tsub">${m(`a_{${id}} = ${T(a)}`)} · segno ${m(`(-1)^{${i + 1}+${j + 1}} = ${s === "+" ? "+1" : "-1"}`)}</span></span>` +
        `<span class="chip" style="background:${color}">${m(`${s}\\, ${PT(a)} \\cdot ${PT(val)} = ${T(contrib)}`)}</span></summary>` +
        `<div class="tbody"><div class="trow">${mat(M, (rr, cc) => (rr === i && cc === j ? "pick" : rr === i || cc === j ? "cut" : ""))}` +
        `<span class="arrow">cancello riga ${i + 1} e colonna ${j + 1}<br>→</span>${mat(sub_)}</div>` +
        (sz > 1
          ? `<div class="nest"><div class="nlabel">Calcolo di ${m(`\\det(M_{${id}})`)} <small>(minore ${sz}×${sz}; gli indici sono riferiti a questa matrice)</small></div>${r.html}</div>`
          : "") +
        `<div class="tres">${m(`\\det(M_{${id}}) = ${T(val)}`)}</div></div></details>`,
    );
  }
  const paint = (t, x) =>
    `${x ? (t.s === "-" ? " - " : " + ") : t.s === "-" ? "-" : ""}\\textcolor{${t.color}}{${PT(t.a)} \\cdot ${PT(t.val)}}`;
  const expr = terms.length ? terms.map((t, x) => paint(t, x)).join("") : "0";
  parts.push(
    `<div class="lsum">Somma dei termini${depth ? "" : " = determinante"}:<br>${m(`\\displaystyle \\det = ${expr} = ${T(total)}`)}</div>`,
  );
  return { det: total, html: parts.join("") };
}

function laplaceSteps(M, pick) {
  const r = laplace(M, pick, true);
  return {
    det: r.det,
    html:
      `<div class="step wide"><div class="stitle">Sviluppo di Laplace</div>` +
      `<div class="ctl"><button type="button" onclick="this.closest('.step').querySelectorAll('details').forEach((d) => (d.open = true))">Espandi tutto</button>` +
      `<button type="button" onclick="this.closest('.step').querySelectorAll('details').forEach((d) => (d.open = false))">Comprimi tutto</button></div>` +
      `<div class="sbody col">${r.html}</div></div>`,
  };
}

/* ---------- Dispatcher ---------- */
const METHODS = {
  gauss: {
    name: "Eliminazione di Gauss",
    ok: (n) => n >= 1,
    run: (M, o = {}) => gauss(M, false, o.partial),
  },
  jordan: {
    name: "Gauss-Jordan",
    ok: (n) => n >= 1,
    run: (M, o = {}) => gauss(M, true, o.partial),
  },
  laplace: {
    name: "Sviluppo di Laplace",
    ok: (n) => n >= 1 && n <= 6,
    run: (M, o = {}) =>
      M.length === 1
        ? {
            det: M[0][0],
            html: step("Matrice 1×1", m("\\det(A) = a_{11}"), mat(M)),
          }
        : laplaceSteps(M, o.pick),
  },
  sarrus: { name: "Regola di Sarrus", ok: (n) => n === 3, run: sarrus },
  det2: { name: "Formula 2×2", ok: (n) => n === 2, run: det2 },
};

if (typeof module !== "undefined") module.exports = { Fr, METHODS };
