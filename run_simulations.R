library("optparse")
suppressWarnings(source("measurement_error_check.R"))

option_list <- list(
  make_option(c("-s", "--simulations"), type = "character", default = NULL, 
              help = "indices of simulations to run, comma separated list of integers and ranges of the n:m, max. range index is 120", metavar = "character"),
  make_option(c("-l", "--label"), type = "character", default = NULL, 
              help = "label for the simulations", metavar = "character"),
  make_option(c("-f", "--filter"), type = "character", default = NULL, 
              help = "Comma separated list of conditions on simulations parameters of form simu_param=value", metavar = "character"),
  make_option(c("-n", "--n"), type = "integer", default = 1, 
              help = "Number of simulations to run", metavar = "integer"),
  make_option(c("-o", "--outdir"), type = "character", default = "data/simulations", 
              help = "output director name [default= %default]", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

parse_simulation_ids <- function(simulations){
  if(is.null(simulations)){
    return(NULL)
  }
  list_simu <- suppressWarnings(trimws(strsplit(simulations, ",")[[1]]))
  ids <- c()
  map(list_simu, function(el){
    range_simu <- suppressWarnings(as.integer(strsplit(el, ":")[[1]]))
    #print(range_simu)
    if(any(is.na(range_simu))){
      stop(sprintf("Invalid range given: '%s'", simulations))  
    }
    if(length(range_simu) != 2){
      range_simu[2] <- range_simu[1]
    }
    ids <<- c(ids, range_simu[1]:range_simu[2])    
  })
  
  #browser() 
  ids <- sort(unique(ids))
  ids[ids > 0 & ids <= 120]
}

parse_filter <- function(filter_def){
  vars <- c("error_level", "size", "with_na", "method")
  var_types <- c("error_level"="character", "size"="integer", "with_na"="logical", "method"="character")
  list_filter <- suppressWarnings(trimws(strsplit(filter_def, ",")[[1]]))
  simu_defs <- get_simu_def()
  filters <- tibble()
  
  map_dfr(list_filter, function(el){
    condition <- suppressWarnings(trimws(strsplit(el, "=")[[1]]))

    if(any(is.na(condition)) || length(condition) != 2){
      stop(sprintf("Invalid filter given: '%s'", filter_def))  
    }
    if(!(condition[1] %in% vars)){
      stop(sprintf("Invalid variable '%s' in filter", condition[1]))  
    }
    allowed_values <- simu_defs %>% pull(!!sym(condition[1])) %>% unique()    

    if(!(condition[2] %in% allowed_values)){
      stop(sprintf("Invalid value '%s' in filter", condition[2]))  
    }
    tibble(var = condition[1], value = condition[2], type = var_types[condition[1]])
  }) 
}

apply_filter <- function(simu_defs, filter_df){
  for(i in 1:nrow(filter_df)){
    if(filter_df[i,]$type == "character"){
      condition <- sprintf("%s=='%s'", filter_df[i,]$var, filter_df[i,]$value)
    }
    else{
      condition <- sprintf("%s==%s", filter_df[i,]$var, filter_df[i,]$value)
      
    }
    simu_defs <- simu_defs %>% filter(eval(str2expression(condition)))
    if(nrow(simu_defs) == 0){
      stop("Filter too strong")
    }
  }
  simu_defs
}

run_simulations <- function(opt){
  simu_ids <- parse_simulation_ids(opt$simulations)
  
  if(!file.exists(opt$outdir)){
    stop(sprintf("Output directory does not exist: '%s'", opt$outdir))  
  }
  if(!is.null(opt$filter)){
    filter_df <- parse_filter(opt$filter)
    simu_defs <- apply_filter(get_simu_def(), filter_df)
  }
  else{
    simu_defs <- get_simu_def() 
  }
  if(!is.null(simu_ids)){
    simu_defs <- simu_defs %>% filter(id %in% simu_ids)
    
  }
  if(nrow(simu_defs) == 0){
    stop("No simulations to execute")
  }
  label <- opt$label
  if(is.null(label)){
    label <- sprintf("simu_%d_%d", min(simu_defs$id), max(simu_defs$id))
  }
  
  setup_workspace_me(F)
  messagef("Running %d simulations range %s-%s",  nrow(simu_defs), min(simu_defs$id), max(simu_defs$id))
  test_all_simulations(master_cross, n_simul = opt$n, simu_def = simu_defs, label = label, out_dir = opt$outdir)
}
#screen -S simu1 
#Rscript --vanilla run_simulations.R  -n 30 -f "size = 250" -l simu_250 -o simulations
#screen -S simu2 
#Rscript --vanilla run_simulations.R  -n 30 -f "size = 500" -l simu_500 -o simulations
#screen -S simu1  
#Rscript --vanilla run_simulations.R  -n 30 -f "size = 1000" -l simu_1000 -o simulations
#screen -S simu2
#Rscript --vanilla run_simulations.R  -n 30 -f "size = 2000" -l simu_2000 -o simulations
#screen -S simu3
#Rscript --vanilla run_simulations.R  -n 30 -f "size = 4000" -l simu_4000 -o simulations

run_simulations(opt)