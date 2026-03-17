document.addEventListener("DOMContentLoaded", () => {
  // RTD theme uses checkbox inputs to manage expansion; force them open.
  document.querySelectorAll("input.toctree-checkbox").forEach((cb) => {
    cb.checked = true;
  });

  // Fallback for any items using expander buttons/spans.
  document.querySelectorAll(".wy-menu-vertical .toctree-expand").forEach((expander) => {
    const li = expander.closest("li");
    if (li) {
      li.classList.add("current");
    }
  });
});
