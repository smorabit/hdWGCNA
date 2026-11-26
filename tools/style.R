# hdWGCNA code style for styler (see STYLE.md).
#
# tidyverse_style() with two changes:
#   - indent by 2 (same as tidyverse)
#   - strings use single quotes: styler's default `fix_quotes` (which converts to double
#     quotes) is replaced by a transformer converting "..." to '...' whenever that is
#     safe (no ' inside the string, no escaped \" inside, not a raw string).
#
# Usage (from the package root):
#   pixi run style        # restyle R/, tests/, data-raw/, bench/
#   or in R: source("tools/style.R"); styler::style_pkg(style = hdwgcna_style)

fix_quotes_single <- function(pd_flat) {
  is_str <- pd_flat$token == 'STR_CONST'
  if (!any(is_str)) {
    return(pd_flat)
  }
  txt <- pd_flat$text[is_str]
  convertible <- startsWith(txt, '"') &
    !grepl("'", txt, fixed = TRUE) &
    !grepl('\\"', txt, fixed = TRUE)
  txt[convertible] <- paste0("'", substr(txt[convertible], 2, nchar(txt[convertible]) - 1), "'")
  pd_flat$text[is_str] <- txt
  pd_flat
}

hdwgcna_style <- function(...) {
  s <- styler::tidyverse_style(indent_by = 2, ...)
  s$token$fix_quotes <- fix_quotes_single
  s$style_guide_name <- 'hdWGCNA::hdwgcna_style@https://github.com/smorabit/hdWGCNA/blob/dev/STYLE.md'
  s
}

style_hdwgcna <- function(paths = c('R', 'tests', 'data-raw', 'bench'), dry = 'off') {
  paths <- paths[dir.exists(paths)]
  out <- lapply(paths, function(p) {
    styler::style_dir(p, transformers = hdwgcna_style(), dry = dry)
  })
  invisible(do.call(rbind, out))
}
