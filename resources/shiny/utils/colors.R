# COLOR PALETTES AND COLOR UTILITY FUNCTIONS -----------------------------------

nice_cols = c(
  "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72", "#B17BA6",
  "#FF7F00", "#FDB462", "#E7298A", "#E78AC3", "#33A02C", "#B2DF8A",
  "#55A1B1", "#8DD3C7", "#A6761D", "#E6AB02", "#7570B3", "#BEAED4",
  "#666666", "#999999", "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3",
  "#808000", "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00",
  "#F8766D", "#E68613", "#7CAE00", "#FF61CC", "#ABA300", "#00C19A",
  "#FEA800", "#02756D", "#FEB479", "#B449BF", "#FE3781", "#52E8F3"
)

palette_map = list(
  "viridis"      = list(option = "viridis", direction =  1, cols = viridis_pal(option = "viridis")(10)),
  "viridis_rev"  = list(option = "viridis", direction = -1, cols = viridis_pal(option = 'viridis', direction = -1)(10)),
  "magma"        = list(option = "magma",   direction =  1, cols = viridis_pal(option = "magma")(10)),
  "magma_rev"    = list(option = "magma",   direction = -1, cols = viridis_pal(option = 'magma', direction = -1)(10)),
  "greytopurple" = list(cols = colorRampPalette(c("#E5E7E9", "#8520B0"))(10)),
  "greytored"    = list(cols = colorRampPalette(c("#E5E7E9", "#EF8611", "#9531C1", "#491C5D"))(10))
)

# Build a colorRamp2 function for use in ComplexHeatmap.
# range controls how min/max are determined:
#   'has_zero'  – one bound is zero; all values must share the same sign
#   'symmetric' – symmetric around zero (|min| == |max|)
#   'natural'   – use the actual data range
cols_fn <- function(mat, res, pal, pal_dir = 1, range = 'natural') {
  stopifnot(pal_dir %in% c(-1, 1))
  range    = match.arg(range, choices = c('has_zero', 'symmetric', 'natural'))
  mat      = mat[!is.na(mat)]
  stopifnot(length(mat) > 0)

  n_vals   = length(unique(mat))
  if (n_vals == 1) {
    pal_cols    = .get_pal_cols(pal, 3)
    cols        = pal_cols[[1]]
    names(cols) = unique(mat)
    return(cols)
  }

  sgn = 1
  if (range == 'has_zero') {
    assert_that( all(mat >= 0) | all(mat <= 0) )
    if (all(mat <= 0)) {
      sgn = -1
      mat = mat * -1
    }
    max_val = mat %>% as.vector %>% max %>% `/`(res) %>% ceiling %>% `*`(res)
    max_val = round(max_val, 3)
    min_val = 0
  } else if (range == 'symmetric') {
    max_val = mat %>% as.vector %>% abs %>% max %>% `/`(res) %>% ceiling %>% `*`(res)
    max_val = round(max_val, 3)
    min_val = -max_val
  } else {
    max_val = mat %>% as.vector %>% max %>% `/`(res) %>% ceiling %>% `*`(res)
    min_val = mat %>% as.vector %>% min %>% `/`(res) %>% floor  %>% `*`(res)
    max_val = round(max_val, 3)
    min_val = round(min_val, 3)
  }

  seq_vals = seq(min_val, max_val, res)
  n_cols   = length(seq_vals)
  pal_cols = .get_pal_cols(pal, n_cols)
  if (length(pal_cols) < n_cols) {
    n_cols   = length(pal_cols)
    seq_vals = seq(min_val, max_val, length.out = n_cols)
  }

  if (pal_dir == -1)
    pal_cols = rev(pal_cols)

  if (sgn == 1)
    cols = colorRamp2(seq_vals, pal_cols)
  else
    cols = colorRamp2(-seq_vals, rev(pal_cols))

  return(cols)
}

.get_pal_cols <- function(pal_str = c('viridis', 'magma', 'Blues', 'BuGn',
                                      'BuPu', 'GnBu', 'Greens', 'Greys', 'Oranges', 'OrRd', 'PuBu', 'PuBuGn',
                                      'PuRd', 'Purples', 'RdPu', 'Reds', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd',
                                      'BrBG', 'PiYG', 'PRGn', 'PuOr', 'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn',
                                      'Spectral', 'Accent', 'Dark2', 'Paired', 'Pastel1', 'Pastel2',
                                      'Set1', 'Set2', 'Set3'), n_cols) {
  pal_str = match.arg(pal_str)
  if (pal_str == 'viridis') {
    pal_cols = str_remove(viridis(n_cols), 'FF$')
  } else if (pal_str == 'magma') {
    pal_cols = str_remove(magma(n_cols), 'FF$')
  } else {
    suppressWarnings({ pal_cols = brewer.pal(n_cols, pal_str) })
  }
  return(pal_cols)
}

# ---- Named palette support --------------------------------------------------

# ggsci palette names supported as flat names
.ggsci_palettes <- c(
  "npg", "aaas", "nejm", "lancet", "jama", "jco", "ucscgb",
  "d3", "d3_10", "d3_20", "d3_20b", "d3_20c",
  "igv", "cosmic", "simpsons", "rickandmorty", "futurama", "tron",
  "startrek", "uchicago", "frontiers", "flatui", "bootstrap"
)

# All valid palette names (built at load time from package metadata)
VALID_PALETTE_NAMES <- c(
  "nice_cols",
  names(MetBrewer::MetPalettes),
  rownames(RColorBrewer::brewer.pal.info),
  .ggsci_palettes
)

# Return n colours from a named palette. Handles overflow via interpolation.
resolve_palette <- function(name, n) {
  stopifnot(is.character(name), length(name) == 1L)
  n <- as.integer(n)
  stopifnot(n > 0L)

  if (name == "nice_cols") {
    return(rep_len(nice_cols, n))
  }

  if (name %in% names(MetBrewer::MetPalettes)) {
    max_n <- length(MetBrewer::MetPalettes[[name]]$colors)
    if (n <= max_n) {
      return(MetBrewer::met.brewer(name, n, type = "discrete"))
    } else {
      return(MetBrewer::met.brewer(name, n, type = "continuous"))
    }
  }

  if (name %in% rownames(RColorBrewer::brewer.pal.info)) {
    max_n     <- RColorBrewer::brewer.pal.info[name, "maxcolors"]
    base_cols <- suppressWarnings(RColorBrewer::brewer.pal(min(n, max_n), name))
    if (n > max_n) return(colorRampPalette(base_cols)(n))
    return(base_cols[seq_len(n)])
  }

  if (name %in% .ggsci_palettes) {
    cols <- .ggsci_cols(name)
    if (n <= length(cols)) return(cols[seq_len(n)])
    return(colorRampPalette(cols)(n))
  }

  stop("Unknown palette '", name, "'. See VALID_PALETTE_NAMES for valid options.")
}

# Internal: retrieve the full discrete colour vector for a ggsci palette name.
.ggsci_cols <- function(name) {
  switch(name,
    "d3"     = ,
    "d3_20"  = ggsci::pal_d3("category20")(20),
    "d3_10"  = ggsci::pal_d3("category10")(10),
    "d3_20b" = ggsci::pal_d3("category20b")(20),
    "d3_20c" = ggsci::pal_d3("category20c")(20),
    "npg"          = ggsci::pal_npg()(10),
    "aaas"         = ggsci::pal_aaas()(10),
    "nejm"         = ggsci::pal_nejm()(8),
    "lancet"       = ggsci::pal_lancet()(9),
    "jama"         = ggsci::pal_jama()(7),
    "jco"          = ggsci::pal_jco()(10),
    "ucscgb"       = ggsci::pal_ucscgb()(26),
    "igv"          = ggsci::pal_igv()(51),
    "cosmic"       = ggsci::pal_cosmic()(26),
    "simpsons"     = ggsci::pal_simpsons()(16),
    "rickandmorty" = ggsci::pal_rickandmorty()(12),
    "futurama"     = ggsci::pal_futurama()(12),
    "tron"         = ggsci::pal_tron()(7),
    "startrek"     = ggsci::pal_startrek()(7),
    "uchicago"     = ggsci::pal_uchicago()(9),
    "frontiers"    = ggsci::pal_frontiers()(10),
    "flatui"       = ggsci::pal_flatui()(9),
    "bootstrap"    = ggsci::pal_bootstrap()(5),
    stop("Unknown ggsci palette: ", name)
  )
}


# Returns "black" or "white" for readable text overlay on a hex background colour.
is_light_or_dark <- function(hex_color) {
  hcl       = as(hex2RGB(hex_color), "polarLUV")
  luminance = hcl[1]
  if (luminance >= 128) "black" else "white"
}
