# hdWGCNA code style

Agreed 2026-11-26. Based on the existing code in `R/` and an interview with Sam.

In short: **tidyverse spacing and layout, single quotes, Seurat-style naming.**
The style is enforced with `styler` (custom style in `tools/style.R`) and `lintr`
(`.lintr`, checked in CI):

```sh
pixi run style   # restyle R/, tests/, data-raw/, bench/
pixi run lint    # report style issues
```

## Layout

- Indent with **2 spaces**, no tabs.
- Lines ≤ 100 characters.
- No trailing whitespace.
- One blank line between logical blocks inside a function, each block introduced by a comment.

## Spacing and braces (tidyverse)

- Space after `if`, `for`, `while` and before `{`: `if (is.null(x)) {`, `for (cur_mod in mods) {`.
- `} else {` and `} else if (...) {` on the closing-brace line.
- Braced bodies go on their own lines, including short guards:

  ```r
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  ```

- Spaces around `=` in arguments: `ScaleData(x, features = genes, verbose = FALSE)`.
- Spaces around `<-`, comparison, arithmetic and logical operators; after commas; none
  inside parentheses.

## Assignment, values, strings

- Assign with `<-`, never `=` or `->`.
- `TRUE` / `FALSE`, never `T` / `F`.
- **Single quotes** for strings: `'grey'`, `'MEs'`. Use double quotes only when the string
  contains a single quote: `"it's"`.
- Use `seq_len()` / `seq_along()` instead of `1:length(x)`, so empty inputs don't break the loop.

## Naming

- **Exported functions:** PascalCase, starting with a verb where possible, in the style of Seurat:
  `ConstructNetwork`, `ModuleEigengenes`, `PlotDMEsVolcano`, `GetModules` / `SetModules`.
- **Internal helpers:** snake_case: `check_wgcna_name()`, `make_pseudobulk_metadata()`.
  This separates them from the exported API at a glance. Rename existing internal helpers
  when touched. (If an internal PascalCase function turns out to be used by people, export
  it rather than rename it.)
- **Arguments:**
  - snake_case for hdWGCNA's own arguments: `wgcna_name`, `soft_power`, `tom_name`,
    `use_metacells`.
  - dot.case only where an argument mirrors a Seurat argument users already know:
    `group.by`, `group.by.vars`, `ident.group`, `vars.to.regress`, `logfc.threshold`.
  - WGCNA pass-through arguments keep WGCNA's camelCase: `networkType`, `minModuleSize`,
    `deepSplit`.
- **Standard argument names** (always the same name for the same thing):
  - `seurat_obj` (first argument), `wgcna_name` (default `NULL` = active experiment)
  - `group.by`, `group_name`, `assay`, `layer`, `verbose`
- **Local variables:** snake_case. Recurring names:
  - `cur_*` for the current element of a loop (`cur_mod`, `cur_genes`, `cur_df`, `cur_seurat`)
  - `plot_df` for data going into a plot, `p` for a ggplot object, `plot_list` for several plots
  - `modules` / `mods` (module table / module names), `MEs`, `kMEs`, `datExpr`, `mod_colors`

## Function structure

```r
ExampleFunction <- function(
  seurat_obj,
  group.by = NULL,
  min_cells = 10,
  wgcna_name = NULL,
  verbose = TRUE,
  ...
) {
  # set the active wgcna experiment if none is given
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  CheckWGCNAName(seurat_obj, wgcna_name)

  # check inputs
  if (!(group.by %in% colnames(seurat_obj@meta.data))) {
    stop(paste0('group.by column ', group.by, ' not found in seurat_obj@meta.data.'))
  }

  # get the modules
  modules <- GetModules(seurat_obj, wgcna_name)

  # ...

  seurat_obj
}
```

- Multi-line signature: `seurat_obj` first, then one argument per line, `...` last.
- Always begin with the `wgcna_name` default + `CheckWGCNAName()` when the function reads
  hdWGCNA data.
- Read and write hdWGCNA data through the getters and setters (`GetModules`, `SetMEs`, ...),
  not `seurat_obj@misc[[...]]` directly.
- **Return the last expression** (`seurat_obj`, `p`); use `return()` only for early exits.
  Functions that modify the object return `seurat_obj`.

## Comments

- Lowercase, short, explaining what the next block does, in the same voice as now:
  `# remove grey module`, `# what features are we using?`, `# transpose MEs so cells are columns`.
- Delete commented-out code; git keeps the history.
- No `# New parameter ...` / `# (Probably fine ...)` notes in finished code. Put them
  in the roxygen docs or the commit message instead.

## Messages and errors

- `stop()` for invalid input, with a message that names the argument and says how to fix it.
- `warning()` for recoverable problems; `message()` for progress, always behind `verbose`.
  No `print()` / `cat()` for user feedback.

## Pipes and tidyverse

- `%>%` (magrittr), not the native `|>`, to keep one pipe style throughout the package.
- dplyr verbs are fine; qualify them (`dplyr::filter`) or import them with `@importFrom`.
- No deprecated tidyverse APIs: `aes_string`, `rename_`, `mutate_at`, ggplot2 `size` for
  lines (use `linewidth`).

## ggplot code

```r
# plot basics
p <- plot_df %>%
  ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj), fill = module))

# optional layers added conditionally
if (show_cutoff) {
  p <- p + geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey75')
}

p <- p +
  geom_point(size = point_size, pch = 21, color = 'black') +
  scale_fill_manual(values = mod_colors) +
  theme(
    panel.border = element_rect(color = 'black', fill = NA, linewidth = 1),
    panel.grid.major = element_blank()
  )

p
```

- Build `p` incrementally, `+` at the end of the line, one layer per line.
- Plot functions return the plot (or a list / patchwork); they don't print or write files.

## Namespaces and dependencies

- `library(hdWGCNA)` attaches **Seurat only** (`Depends`). Everything else is in `Imports`
  or `Suggests`.
- Functions from `Imports` packages: either `@importFrom pkg fn` or `pkg::fn()`. No
  whole-package `@import`.
- `Suggests` packages (optional features): always `pkg::fn()`, guarded at the top of the
  function:

  ```r
  if (!requireNamespace('DESeq2', quietly = TRUE)) {
    stop("Package 'DESeq2' is required for this function. Install it with BiocManager::install('DESeq2').")
  }
  ```

- Never `library()`, `require()`, or install packages inside package code.

## Documentation (roxygen2)

Depth scales with how much users interact with the function:

- **Main user-facing functions** (e.g. `SetupForWGCNA`, `MetacellsByGroups`, `SetDatExpr`,
  `TestSoftPowers`, `ConstructNetwork`, `ModuleEigengenes`, `ModuleConnectivity`, `FindDMEs`,
  `ModulePreservation`, plotting functions): full documentation:
  - title + `@description` paragraph explaining what it does and when to use it
  - `@param` for every argument with its type, meaning, and default behaviour
  - `@details` for the method, prerequisites (which functions must run first), and caveats
  - `@return` describing exactly what is added to the object / returned
  - `@examples` running on `hdw_sim` (fast) or inside `\donttest{}`
  - `@seealso` linking related functions
- **Small functions** (getters/setters, simple utilities): short docs: a title, one-line
  `@param`s, `@return`.
- **Internal helpers:** `@noRd` with a short comment, or plain comments.
- Every exported function has `@return` and `@export`.
- `@keywords`: one term each (the `@keywords scRNA-seq, plot,` form breaks the Rd index),
  or drop them.

## Tests

- One file per source area: `tests/testthat/test-<Function or topic>.R`.
- Use the shared fixture built on `hdw_sim`; write files only to `tempdir()`.
- Test names describe the behaviour: `test_that('ModuleEigengenes returns one ME per module', ...)`.

## Applying the style

- Phase 1 includes one formatting-only commit (`pixi run style`) so later diffs stay
  readable. Use `git blame --ignore-rev <sha>` (or a `.git-blame-ignore-revs` file) to skip it.
- Naming changes (helper renames, argument harmonisation) are **not** done by the formatter;
  they happen by hand during the cleanup, with deprecation shims for exported arguments.
- The lint workflow reports only until the cleanup is done, then fails on any lint.
