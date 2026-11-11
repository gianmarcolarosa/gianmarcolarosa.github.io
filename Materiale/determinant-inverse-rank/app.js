"use strict";
// UI: input grid, method choice, step rendering.

const $ = (id) => document.getElementById(id);
let vals = [
  ["2", "1", "3"],
  ["0", "-1", "4"],
  ["1", "2", "5"],
];
let method = "gauss";

const dims = () => [clampInt($("rows").value), clampInt($("cols").value)];
const clampInt = (v) => Math.min(8, Math.max(1, parseInt(v, 10) || 1));

function buildGrid() {
  const [r, c] = dims();
  vals = Array.from({ length: r }, (_, i) =>
    Array.from({ length: c }, (_, j) => (vals[i] && vals[i][j]) || "0"),
  );
  $("grid").innerHTML = `<div class="matw big"><table class="mat">${vals
    .map(
      (row, i) =>
        `<tr>${row.map((v, j) => `<td><input data-i="${i}" data-j="${j}" value="${v}" inputmode="decimal"></td>`).join("")}</tr>`,
    )
    .join("")}</table></div>`;
  $("grid")
    .querySelectorAll("input")
    .forEach((el) =>
      el.addEventListener("input", () => {
        vals[+el.dataset.i][+el.dataset.j] = el.value;
        el.classList.remove("bad");
      }),
    );
  refreshMethods();
}

function refreshMethods() {
  const [r, c] = dims();
  const square = r === c;
  $("msg").textContent = square
    ? ""
    : `Matrice ${r}×${c}: determinante e inversa sono definiti solo per matrici quadrate (righe = colonne); il rango no.`;
  $("go").disabled = !square;
  $("methods").innerHTML = Object.entries(METHODS)
    .map(([id, m]) => {
      const ok = square && m.ok(r);
      const why = !square
        ? ""
        : id === "sarrus"
          ? "solo 3×3"
          : id === "det2"
            ? "solo 2×2"
            : id === "laplace"
              ? "max 6×6"
              : "";
      return (
        `<label class="method ${ok ? "" : "off"}"><input type="radio" name="m" value="${id}" ${ok ? "" : "disabled"} ${id === method ? "checked" : ""}>` +
        `${m.name}${why ? ` <small>(${why})</small>` : ""}</label>`
      );
    })
    .join("");
  if (square && !METHODS[method].ok(r)) {
    method = "gauss";
    const first = $("methods").querySelector('input[value="gauss"]');
    if (first) first.checked = true;
  }
  $("methods")
    .querySelectorAll("input")
    .forEach((el) =>
      el.addEventListener("change", () => {
        method = el.value;
        showLap();
      }),
    );
  showLap();
}

function showLap() {
  const [r] = dims();
  $("lapopt").hidden = method !== "laplace";
  $("pivopt").hidden = method !== "gauss" && method !== "jordan";
  const t = $("laptype").value,
    cur = $("lapidx").value;
  $("lapidx").innerHTML = Array.from(
    { length: r },
    (_, i) =>
      `<option value="${i}">${t === "row" ? "riga" : "colonna"} ${i + 1}</option>`,
  ).join("");
  if (cur && +cur < r) $("lapidx").value = cur;
  $("lapidx").disabled = $("laptype").disabled = $("lapauto").checked;
}

function readMatrix() {
  let ok = true;
  const M = vals.map((row, i) =>
    row.map((v, j) => {
      const f = Fr.parse(v);
      if (!f) {
        ok = false;
        const el = $("grid").querySelector(
          `input[data-i="${i}"][data-j="${j}"]`,
        );
        if (el) el.classList.add("bad");
      }
      return f;
    }),
  );
  return ok ? M : null;
}

function calculate() {
  const M = readMatrix();
  if (!M) {
    $("msg").textContent =
      "Alcuni valori non sono numeri validi (evidenziati in rosso).";
    return;
  }
  $("msg").textContent = "";
  const n = M.length,
    m = METHODS[method];
  const pick =
    method === "laplace" && !$("lapauto").checked
      ? { t: $("laptype").value, i: +$("lapidx").value }
      : null;
  const res = m.run(M, { pick, partial: $("pivot").value === "partial" });
  const approx = res.det.d !== 1n ? ` <small>≈ ${res.det.dec()}</small>` : "";
  $("result").innerHTML =
    `<section class="card"><h2>3 · ${m.name} — matrice ${n}×${n}</h2><div class="steps">${res.html}</div>` +
    `<div class="final">det(A) = ${res.det.html()}${approx}</div></section>`;
  // ponytail: needs internet for KaTeX; formulas show as raw LaTeX offline
  if (window.renderMathInElement)
    renderMathInElement($("result"), {
      delimiters: [{ left: "\\(", right: "\\)", display: false }],
    });
  $("result").scrollIntoView({ behavior: "smooth", block: "start" });
}

$("rows").addEventListener("change", buildGrid);
$("cols").addEventListener("change", buildGrid);
$("go").addEventListener("click", calculate);
$("laptype").addEventListener("change", showLap);
$("lapauto").addEventListener("change", showLap);
$("example").addEventListener("click", () => {
  const [r, c] = dims();
  vals = Array.from({ length: r }, (_, i) =>
    Array.from({ length: c }, (_, j) =>
      String(((i * 3 + j * 5 + i * j) % 7) - 2),
    ),
  );
  buildGrid();
});
$("random").addEventListener("click", () => {
  const [r, c] = dims();
  vals = Array.from({ length: r }, () =>
    Array.from({ length: c }, () => String(Math.floor(Math.random() * 11) - 5)),
  );
  buildGrid();
});
$("clear").addEventListener("click", () => {
  vals = [];
  buildGrid();
});

buildGrid();
