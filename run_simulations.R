library("optparse")
suppressWarnings(source("measurement_error_check.R"))

option_list <- list(
  make_option(c("-s", "--simulations"), type = "character", default = NULL, 
              help = "indices of simulations to run, form n:m, max. range 1:108", metavar = "character"),
  make_option(c("-n", "--n"), type = "integer", default = 1, 
              help = "Number of simulations to run", metavar = "integer"),
  make_option(c("-o", "--outdir"), type = "character", default = "data/simulations", 
              help = "output director name [default= %default]", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

run_simulations <- function(opt){
  if(is.null(opt$simulations)){
    stop("No simulations provided")
  }   
  range_simu <- suppressWarnings(as.integer(strsplit(opt$simulations, ":")[[1]]))
  if(any(is.na(range_simu)) || length(range_simu) != 2){
    stop(sprintf("Invalid range given: '%s'", opt$simulations))  
  }  
  if(!file.exists(opt$outdir)){
    stop(sprintf("Output directory does not exist: '%s'", opt$outdir))  
    
  }
  simu_defs <- get_simu_def()[range_simu[1]:range_simu[2],]
  label <- sprintf("simu_%d_%d", range_simu[1], range_simu[2])
  setup_workspace_me(F)
  test_all_simulations(master_cross, n_simul = opt$n, simu_def = simu_defs, label = label, out_dir = opt$outdir)
}

run_simulations(opt)