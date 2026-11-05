"use strict";
// Theme choice (auto / light / dark), remembered in localStorage; "auto" follows the system.
(function () {
  const KEY = "determinant-inverse-theme";
  const apply = (v) => {
    if (v === "light" || v === "dark") document.documentElement.dataset.theme = v;
    else delete document.documentElement.dataset.theme;
  };
  let saved = "auto";
  try { saved = localStorage.getItem(KEY) || "auto"; } catch (e) { /* file:// without storage */ }
  apply(saved);
  document.addEventListener("DOMContentLoaded", () => {
    const sel = document.getElementById("theme");
    sel.value = saved;
    sel.addEventListener("change", () => {
      apply(sel.value);
      try { localStorage.setItem(KEY, sel.value); } catch (e) { /* ignore */ }
    });
  });
})();
