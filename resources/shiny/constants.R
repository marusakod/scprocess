# CONSTANTS --------------------------------------------------------------------

# Font sizes
FONT_TITLE = 18L
FONT_AXIS  = 16L
FONT_TEXT  = 14L
FONT_SMALL = 12L

# Log-scale axis breaks and labels (for pseudobulk count plots)
LOG_BRKS = c(0, 10, 20, 50, 1e2, 2e2, 5e2, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4) %>% `+`(10) %>% log
LOG_LABS = c('0', '10', '20', '50', '100', '200', '500', '1k', '2k', '5k', '10k', '20k', '50k')

# Log10-scale axis breaks and labels (for composition / prevalence plots)
PCT_BRKS = c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1) %>% log10
PCT_LABS = c('0.1%', '0.2%', '0.5%', '1%', '2%', '5%', '10%', '20%', '50%', '100%')
