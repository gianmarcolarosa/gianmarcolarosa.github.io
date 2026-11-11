"use strict";
// Rank tab: runs the chosen method of rank.js on the matrix of the shared grid and shows the steps.

const rankStepper = INV.makeStepper("st_rank", 1400);
const RANK_METHODS = {
  gauss: ["Eliminazione di Gauss", RANK.gaussSteps],
  gj: ["Gauss–Jordan", RANK.gaussJordanSteps],
  minors: ["Minori orlati", RANK.minorsSteps],
};

function calculateRank() {
  const M = readMatrix();
  if (!M) {
    $("msg").textContent = "Alcuni valori non sono numeri validi (evidenziati in rosso).";
    return;
  }
  $("msg").textContent = "";
  const [label, run] = RANK_METHODS[document.querySelector('input[name="rkm"]:checked').value];
  const A = toInv(M);
  const res = run(A);
  rankStepper.set(res.steps);
  $("rankTitle").textContent = `${label}: rango ${res.rank} (matrice ${A.length}×${A[0].length})`;
  $("rankresult").hidden = false;
  $("rankresult").scrollIntoView({ behavior: "smooth", block: "start" });
}

$("goRank").addEventListener("click", calculateRank);
document.querySelectorAll('input[name="rkm"]').forEach((el) =>
  el.addEventListener("change", () => {
    if (!$("rankresult").hidden) calculateRank();
  }),
);
