"use strict";
// Rank of a matrix, three methods with exact fractions: Gauss (row echelon form), Gauss-Jordan (reduced form) and
// bordering minors (Kronecker). Built on the exact arithmetic and the step format of inverse.js (INV).
window.RANK = (function () {
  const { ZERO, mCopy, mat, T, opSwap, opElim, det, gjEngine } = INV;
  const cells = (rows, cols) => rows.flatMap((i) => cols.map((j) => [i, j]));
  const list = (idx) => idx.map((i) => i + 1).join(", ");

  /** The last step of every method: the rank, the nullity and what they say about the matrix. */
  function conclusion(A, r, pivots, how) {
    const m = A.length, n = A[0].length;
    const cols = pivots
      ? ` Le colonne ${list(pivots)} della matrice <b>originale</b> (quelle dei pivot) sono linearmente indipendenti e generano lo spazio delle colonne: ne sono una base.`
      : "";
    const full = r === Math.min(m, n) ? ` Il rango è massimo (${T(`r = \\min(${m},${n})`)}).` : "";
    const inv = m === n ? ` La matrice è quadrata: ${r === n ? `${T("\\det A \\ne 0")}, quindi è <b>invertibile</b>` : `${T("\\det A = 0")}, quindi <b>non è invertibile</b>`}.` : "";
    return {
      text: `${how} Il <b>rango</b> è ${T(`\\operatorname{rg}(A) = ${r}`)}. Per il teorema rango–nullità la dimensione del nucleo è ${T(`n - r = ${n} - ${r} = ${n - r}`)}.${cols}${full}${inv}`,
      mats: [mat("A", A)],
      vars: { rango: `${r}`, "\\text{nullità}": `${n - r}` },
      tone: "ok",
    };
  }

  /** Method 1: Gauss, only below the pivots; the rank is the number of non-zero rows (pivots). */
  function gaussSteps(A) {
    const m = A.length, n = A[0].length, W = mCopy(A), pivots = [];
    const steps = [{ text: `Porto $A$ in <b>forma a scala</b> con operazioni sulle righe (scambio due righe, sommo a una riga un multiplo di un'altra): queste operazioni non cambiano il rango. In forma a scala il rango è il numero di <b>pivot</b>, cioè di righe non nulle.`, mats: [mat("A", mCopy(W))], vars: { rango: "?" } }];
    const snap = (text, hl, vars) => steps.push({ text, mats: [mat("A", mCopy(W), { hl })], vars });
    let r = 0;
    for (let c = 0; c < n && r < m; c++) {
      let p = -1;
      for (let i = r; i < m; i++) if (!W[i][c].zero) { p = i; break; }
      if (p < 0) {
        snap(`Colonna ${c + 1}: dalla riga ${r + 1} in giù sono tutti $0$, non c'è pivot in questa colonna. Passo alla colonna successiva restando sulla riga ${r + 1}.`, { cols: [c] }, { colonna: `${c + 1}` });
        continue;
      }
      if (p !== r) {
        [W[p], W[r]] = [W[r], W[p]];
        snap(`Colonna ${c + 1}: l'elemento in posizione (${r + 1},${c + 1}) è $0$. Scambio ${T(opSwap(r, p))} per avere un pivot, ${T(W[r][c].tex())}.`, { rows: [r, p], strong: [[r, c]] }, { colonna: `${c + 1}`, pivot: W[r][c].tex(), operazione: opSwap(r, p) });
      }
      for (let i = r + 1; i < m; i++) {
        if (W[i][c].zero) continue;
        const k = W[i][c].div(W[r][c]);
        W[i] = W[i].map((x, j) => x.sub(k.mul(W[r][j])));
        snap(`Azzero l'elemento della riga ${i + 1} sotto il pivot ${T(W[r][c].tex())}: ${T(opElim(i, r, k))}.`, { rows: [i], cols: [c], strong: [[r, c], [i, c]] }, { colonna: `${c + 1}`, pivot: W[r][c].tex(), operazione: opElim(i, r, k) });
      }
      pivots.push(c);
      r++;
    }
    const last = conclusion(A, pivots.length, pivots, `La matrice è in forma a scala: ci sono ${pivots.length} pivot (colonne ${list(pivots) || "nessuna"}).`);
    last.mats = [mat("A", mCopy(A)), mat("\\text{forma a scala}", mCopy(W), { hl: { strong: pivots.map((c, i) => [i, c]) } })];
    steps.push(last);
    return { steps, rank: pivots.length };
  }

  /** Method 2: Gauss-Jordan to the reduced form (pivots equal to 1, zeros above and below). */
  function gaussJordanSteps(A) {
    const n = A[0].length;
    const e = gjEngine(A, { nc: n, partial: false, name: "A" });
    const steps = [{ text: `Porto $A$ in <b>forma ridotta</b> (Gauss–Jordan): ogni pivot diventa $1$ e sopra e sotto di esso ci sono solo $0$. Il rango è il numero di pivot, e le colonne dei pivot dicono quali colonne di $A$ sono indipendenti.`, mats: [mat("A", mCopy(A))], vars: { rango: "?" } }, ...e.steps];
    const last = conclusion(A, e.pivots.length, e.pivots, `La forma ridotta ha ${e.pivots.length} pivot (colonne ${list(e.pivots) || "nessuna"}).`);
    last.mats = [mat("A", mCopy(A)), mat("\\text{forma ridotta}", mCopy(e.W), { hl: { strong: e.pivots.map((c, i) => [i, c]) } })];
    steps.push(last);
    return { steps, rank: e.pivots.length };
  }

  /** Method 3: bordering minors (Kronecker): grow a non-zero minor until every minor that borders it is zero. */
  function minorsSteps(A) {
    const m = A.length, n = A[0].length, top = Math.min(m, n);
    const sub = (I, J) => I.map((i) => J.map((j) => A[i][j]));
    const steps = [{ text: `Il rango è l'ordine massimo di un <b>minore non nullo</b> (determinante di una sottomatrice quadrata). Per non provarli tutti si usa il <b>teorema degli orlati</b>: se un minore $M$ di ordine $r$ è non nullo e <b>tutti</b> i minori di ordine $r+1$ che lo contengono (aggiungendo una riga e una colonna) sono nulli, allora il rango è $r$.`, mats: [mat("A", mCopy(A))], vars: { rango: "?" } }];
    let I = null, J = null;
    for (let i = 0; i < m && !I; i++) for (let j = 0; j < n && !I; j++) if (!A[i][j].zero) { I = [i]; J = [j]; }
    if (!I) {
      const last = conclusion(A, 0, null, `Tutti gli elementi sono $0$: non esiste nessun minore non nullo.`);
      steps.push(last);
      return { steps, rank: 0 };
    }
    steps.push({ text: `Cerco un elemento non nullo: ${T(`a_{${I[0] + 1}${J[0] + 1}} = ${A[I[0]][J[0]].tex()} \\ne 0`)} è un minore di ordine $1$ non nullo, quindi ${T("\\operatorname{rg}(A) \\ge 1")}.`, mats: [mat("A", mCopy(A), { hl: { soft: cells(I, J) } })], vars: { rango: "\\ge 1", "\\text{minore}": A[I[0]][J[0]].tex() } });
    while (I.length < top) {
      const k = I.length + 1, tried = [];
      let found = null;
      for (let i = 0; i < m && !found; i++) {
        if (I.includes(i)) continue;
        for (let j = 0; j < n && !found; j++) {
          if (J.includes(j)) continue;
          const RI = [...I, i].sort((a, b) => a - b), CJ = [...J, j].sort((a, b) => a - b), d = det(sub(RI, CJ));
          tried.push(T(`\\det A_{\\{${list(RI)}\\},\\{${list(CJ)}\\}} = ${d.tex()}`));
          if (!d.zero) found = { i, j, RI, CJ, d };
        }
      }
      const head = `Provo i minori di ordine ${k} che orlano il minore attuale (righe ${list(I)}, colonne ${list(J)}): `;
      if (found) {
        steps.push({ text: `${head}${tried.join(" · ")}. L'ultimo è <b>non nullo</b>: il rango è almeno ${k}.`, mats: [mat("A", mCopy(A), { hl: { soft: cells(I, J), strong: [...found.RI.map((i) => [i, found.j]), ...found.CJ.map((j) => [found.i, j])] } }), mat("\\text{nuovo minore}", sub(found.RI, found.CJ))], vars: { rango: `\\ge ${k}`, "\\text{minore}": found.d.tex() } });
        I = found.RI; J = found.CJ;
      } else {
        const last = conclusion(A, I.length, null, `${head}${tried.join(" · ")}. Sono <b>tutti nulli</b>: per il teorema degli orlati il rango è ${I.length}.`);
        last.mats = [mat("A", mCopy(A), { hl: { soft: cells(I, J) } }), mat("\\text{minore non nullo}", sub(I, J))];
        steps.push(last);
        return { steps, rank: I.length };
      }
    }
    const last = conclusion(A, I.length, null, `Il minore ha ordine ${top} $= \\min(m,n)$: non ci sono minori più grandi.`);
    last.mats = [mat("A", mCopy(A), { hl: { soft: cells(I, J) } }), mat("\\text{minore non nullo}", sub(I, J))];
    steps.push(last);
    return { steps, rank: I.length };
  }

  return { gaussSteps, gaussJordanSteps, minorsSteps };
})();
