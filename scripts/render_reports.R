# render_reports.R
suppressPackageStartupMessages({
  library('stringr')
  library('assertthat')
  library('workflowr')
  library('glue')
})

# function that replaces placeholder strings in template Rmds and renders html reports

# commented this out because we don't use lists with rmd files
# render_reports <- function(proj_dir, rmd_ls_concat) {
#   rmd_ls    = rmd_ls_concat %>% str_split(" ") %>% unlist
#   print(rmd_ls)
#   print(proj_dir)
#   setwd(proj_dir)
#   assert_that( all(file.exists(rmd_ls)) )

#   for (rmd_f in rmd_ls)
#     workflowr::wflow_build( files = rmd_f, view = FALSE,
#       verbose = TRUE, delete_cache = TRUE )
# }

render_reports <- function(proj_dir, temp_f, temp_ls, rmd_f){
 
setwd(proj_dir) 
# make Rmd file
message('Creating Rmd file from template ', temp_f)
make_rmd_from_temp(temp_f, temp_ls, rmd_f)

# render html

# workflowr::wflow_build( files = rmd_f, view = FALSE,
#  verbose = TRUE, delete_cache = TRUE)
return(NULL)


}



make_rmd_from_temp <- function(temp_f, temp_ls, rmd_f){
 if(!file.exists(rmd_f)){
  # read remplate file
  temp_str = readLines(temp_f, warn = FALSE)
  
  # convert into one long string
  temp_str = paste(temp_str, collapse = "\n")

  # substitute placeholders
  rmd_str = glue(temp_str, .envir = as.environment(temp_ls), .open = "${", .close = "}")

  # write to rmd file
  writeLines(rmd_str, rmd_f)
 }


}
