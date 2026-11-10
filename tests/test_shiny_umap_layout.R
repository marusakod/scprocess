suppressPackageStartupMessages({
  library(data.table)
  library(magrittr)
  library(ggplot2)
  library(scattermore)
  library(patchwork)
  library(ggrepel)
})
source('resources/shiny/constants.R')
source('resources/shiny/utils/plots_clusters.R')
source('resources/shiny/utils/plots_genes.R')

# Compare a small atlas with a larger one on the same display, then give the
# larger atlas more height. The legend must adapt rather than always use two columns.
labels = sprintf('Population %02d', seq_len(47))
colours = setNames(grDevices::hcl.colors(47, palette = 'Dark 3'), labels)
centroids = data.table(cluster = factor(labels, levels = labels),
                       UMAP_1 = cos(seq(0, 2 * pi, length.out = 47)),
                       UMAP_2 = sin(seq(0, 2 * pi, length.out = 47)),
                       label_num = seq_len(47), text_color = 'black')
small = cluster_umap_legend(centroids[1:8], colours, c(1000, 680) / 72)
large = cluster_umap_legend(centroids, colours, c(1000, 680) / 72)
tall = cluster_umap_legend(centroids, colours, c(1000, 1000) / 72)
stopifnot(small$n_columns == 1, small$font_size == 12,
          large$n_columns == 2, large$font_size < small$font_size,
          large$font_size >= 9, tall$n_columns == 1)

# Check legend order, palette alignment and bounds on the physical device.
check_legend = function(legend, expected, device_size) {
  symbols = legend$grob$children[[2]]
  numbers = legend$grob$children[[3]]
  text = legend$grob$children[[4]]
  stopifnot(identical(numbers$label, as.character(expected$label_num)),
            identical(symbols$gp$fill, unname(colours[as.character(expected$cluster)])),
            identical(gsub('\n', ' ', text$label, fixed = TRUE), as.character(expected$cluster)),
            all(as.numeric(text$x) > 0), all(as.numeric(text$x) < legend$width),
            all(as.numeric(text$y) > 0), all(as.numeric(text$y) < device_size[2]))
}
check_legend(large, centroids, c(1000, 680) / 72)
long_centroids = copy(centroids)
long_centroids[, cluster := paste('Cortical neuronal population with a long annotation', label_num)]
long_colours = setNames(unname(colours), long_centroids$cluster)
long_legend = cluster_umap_legend(long_centroids, long_colours, c(1000, 680) / 72)
stopifnot(any(grepl('\n', long_legend$grob$children[[4]]$label, fixed = TRUE)),
          long_legend$font_size >= 9, long_legend$font_size <= large$font_size)

# Render on devices matching Shiny's resolution, opening each device before
# constructing the plot. Exercise centroid and repel modes on narrow/wide views.
out = Sys.getenv('SHINY_LAYOUT_TEST_OUTPUT', tempdir())
set.seed(7)
cells = centroids[rep(seq_len(.N), each = 12), .(cluster, UMAP_1, UMAP_2)]
cells[, `:=`(UMAP_1 = UMAP_1 + rnorm(.N, sd = 0.04),
             UMAP_2 = UMAP_2 + rnorm(.N, sd = 0.04))]
repel_pos = centroids[, .(cluster, repel_1 = UMAP_1 * 1.1, repel_2 = UMAP_2 * 1.1)]
for (width in c(500, 1000)) {
  layout = cluster_umap_legend(centroids, colours, c(width, 680) / 72)
  check_legend(layout, centroids, c(width, 680) / 72)
  for (mode in c('centroid', 'repel')) {
    png(file.path(out, paste0('scprocess-layout-', width, '-', mode, '.png')),
        width = width, height = 680, res = 72)
    p = plot_cluster_umap(cells, colours, centroids, repel_pos, label_mode = mode)
    stopifnot(p[[1]]$theme$aspect.ratio == 1)
    print(p)
    dev.off()
  }
}
palette_map = list(viridis = list(option = 'viridis', direction = 1))
cells[, gene_exp := log1p(seq_len(.N) %% 10)]
p = plot_gene_umap(cells, 'Test gene', 'viridis')
stopifnot(p$theme$aspect.ratio == 1)
ggsave(file.path(out, 'scprocess-layout-gene.png'), p,
       width = 1000 / 72, height = 680 / 72, dpi = 72, bg = 'white')
message('PASS: adaptive columns/font size, wrapped labels, order/colours/bounds, square panels and both label modes')
