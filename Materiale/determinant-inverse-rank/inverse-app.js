"use strict";
// Inverse tab: cofactor/adjugate and Gauss-Jordan steps (inverse.js) for the matrix typed in the shared grid.

const cofStepper = INV.makeStepper("st_cof", 1400);
const gjStepper = INV.makeStepper("st_gj", 1400);

// det.js and inverse.js have separate Fr classes: copy numerator and denominator
const toInv = (M) => M.map((row) => row.map((f) => INV.fr(f.n, f.d)));

function showTab(id) {
  for (const t of ["det", "inverse", "rank"]) {
    $("tab-" + t).hidden = t !== id;
    $("t-" + t).setAttribute("aria-selected", String(t === id));
  }
  if (id !== "det") INV.typeset($("tab-" + id)); // typeset the static $...$ of the tab
}

const inverseMethod = () => document.querySelector('input[name="invm"]:checked').value;

function calculateInverse() {
  const M = readMatrix();
  if (!M) {
    $("msg").textContent = "Alcuni valori non sono numeri validi (evidenziati in rosso).";
    return;
  }
  if (M.length !== M[0].length) {
    $("msg").textContent = `Matrice ${M.length}×${M[0].length}: l'inversa è definita solo per matrici quadrate.`;
    return;
  }
  $("msg").textContent = "";
  const A = toInv(M);
  const gj = inverseMethod() === "gj";
  $("sec_cof").hidden = gj;
  $("sec_gj").hidden = !gj;
  if (gj) gjStepper.set(INV.gaussJordanSteps(A, $("gjPartial").checked).steps);
  else cofStepper.set(INV.cofactorSteps(A).steps);
  $("invresult").hidden = false;
  $("invresult").scrollIntoView({ behavior: "smooth", block: "start" });
}

$("t-det").addEventListener("click", () => showTab("det"));
$("t-inverse").addEventListener("click", () => showTab("inverse"));
$("t-rank").addEventListener("click", () => showTab("rank"));
$("goInv").addEventListener("click", calculateInverse);
$("gjPartial").addEventListener("change", () => {
  if (!$("invresult").hidden) calculateInverse();
});
document.querySelectorAll('input[name="invm"]').forEach((el) =>
  el.addEventListener("change", () => {
    $("gjopt").hidden = inverseMethod() !== "gj";
    if (!$("invresult").hidden) calculateInverse();
  }),
);
