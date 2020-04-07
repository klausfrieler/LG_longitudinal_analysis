library(tidyverse)

#'Helper functions that sets a certain @variable @value for a certin @p_id
#'
set_val_for_pid <- function(data, p_id, var_name, value){
  #browser()

  #messagef("%s, %s, %s", p_id, var_name, value)
  if(!(p_id %in% data$p_id)){
    stop("Invalid p_id: %s", p_id)
  }
  data[!is.na(data$p_id) & data$p_id == p_id, var_name] <- value
  data
}

#'Gender is assumed to be constant over a lifetime
#'just to simplify analysis, so if a participant indicated more than one gender 
#'in several test iterations, the "majority" gender will be estimated
#'and stored in the variable "patch_gender"
#'
get_majority_gender <- function(data){
  tmp <- data %>% 
    count(gender) %>% 
    filter(!is.na(gender)) %>% 
    arrange(desc(n)) %>% mutate(r = nrow(.) - rank(n) + 1)
  if(nrow(tmp) == 0){
    return(NA)
  }
  if(nrow(tmp) == 1){
    return(tmp[1,]$gender)
  }
  tmp1 <- tmp %>% filter(gender %in% c("Male", "Female"))
  if(nrow(tmp1) == 1){
    return(tmp1[1, ]$gender)
  }
  if(tmp[1,]$n > tmp[2,]$n){
    return(tmp[1,]$gender)
  }
  return(NA)
}
#'Gender is assumed to be constant over a lifetime
#'just to simplify analysis, so if a participant indicated more than one gender 
#'in several test iterations, the "majority" gender will be estimated
#'and stored in the variable "patch_gender"
#'
patch_gender <- function(data, substitute = F){
  ids <- unique(data$p_id)
  ids <- ids[!is.na(ids)]
  pg <- 
    map_dfr(ids, function(x) { 
    #browser()
    gender <- data %>% filter(p_id == x) %>% get_majority_gender()
    tibble(p_id = x, patch_gender = gender)
  })
  if(substitute){
    data$gender <- NULL
    data <- data %>% left_join(pg %>% select(gender = patch_gender), by = "p_id")
  }
  else{
    data <- data %>% left_join(pg, by = "p_id")
  }
  data
}

#'Age in month can be calculated from birth year and month, if not already set
#'
get_age_months <- function(birth_year, birth_month, ref_time){
  if(is.na(birth_year)){
    return(NA)
  }
  if(is.na(birth_month)){
    birth_month <- 7
  }
  date_of_birth <- lubridate::ymd(sprintf("%s-%02d-15", birth_year, birth_month))
  if(is.character(ref_time)){
    ref_time <- lubridate::ymd(ref_time)
  } 
  interval_period = lubridate::interval(date_of_birth, ref_time )
  #browser()
  full_year = interval_period %/% lubridate::years(1)
  remaining_months = interval_period %% lubridate::years(1) %/% months(1)
  #printf('Your age is %d years, %d month', full_year, remaining_months)
  full_year * 12 + remaining_months
}

#' Some earlier GMSI versions used sum scores instead of mean scores
#' 
rescale_GMSI <- function(df) {
  #browser()
  items <- read.csv("data/gms-items.csv", stringsAsFactors = FALSE)
  num_items <- tolower(items$factor_label) %>% gsub(" ", "_", .) %>% table
  num_items[["general"]] <- sum(items$general_factor)
  for (factor in names(num_items)) {
    label <- paste0("GMS.", factor)
    to_update <- which(df[[label]] > 7)
    n <- num_items[[factor]]
    if(length(label) == 0 & length(to_update) == 0){
      browser()
    }
    df[[label]][to_update] <- df[[label]][to_update] / n
    stopifnot(all(df[[label]] >= 1 & df[[label]] <= 7, na.rm = TRUE))
  }
  df
}
#'Year group is the term for year groups (UK) and "Klassenstufe" (DE)
#'
make_year_group <-  Vectorize( 
  function(year_group_uk, year_group_de, age, country ){
    na_uk <- is.na(year_group_uk)
    na_de <- is.na(year_group_de)
    
    if(na_uk && na_de){
      year_group <- year_group_int_from_age(age, country)
      #messagef("YG DE: %s, YG EK:%s, age: %d, country: %s -> %s", year_group_de, year_group_uk, age, country, year_group)
      return(year_group)
    }
    if(!is.na(country)){
      year_group <- ifelse(country == "DE", year_group_de, year_group_uk)
      if(is.na(year_group)){
        year_group <- sum(c(year_group_de, year_group_uk), na.rm = T)
      }
    }
    else{
      year_group <- NA
      #messagef("YG DE: %s, YG EK:%s, age: %d, country: %s -> %s", year_group_de, year_group_uk, age, country, year_group)
    }
    year_group
  }
)

age_group_from_age <- Vectorize( 
  function(age, country){
    if(is.na(age) || age < 10 || age > 18 || is.na(country)){
      return(NA) 
    }
    age_groups <- sprintf("%d_%d", 10:17, 11:18)
    age_groups <- c(age_groups, "17_18")
    age_groups[age - 9]
    
  }
)

year_group_int_from_age <- Vectorize( 
  function(age, country){
    if(is.na(age) ||  is.na(country) || !(country %in% c("DE", "UK"))){
      return(NA) 
    }
    if(country == "DE"){
      year_group <- age - 5
      year_group <- max(year_group, 5)
    }  
    if(country == "UK"){
      year_group <- age - 4
      year_group <- max(year_group, 7)
    }
    year_group <- min(year_group, 13)
    year_group
  }
)
#' Age group is a transnational hybrid of year group and age
#' to reflect the differente age structure of UK and German students,
#' as the UK students enter school a year earlier than the Germans
#' So, a 5th grade in UK is typically 11 years, whereas in German they are 12 years old
#' Age group combines UK n-grader and German n+1-grader in the same grop.
#'  
make_age_group <- Vectorize(
  function(year_group_uk, year_group_de, age, country){
    #messagef("YG UK: %s, YG DE: %s, age: %s, country: %s", year_group_uk, year_group_de, age, country)
    na_uk <- is.na(year_group_uk)
    na_de <- is.na(year_group_de)
    if(is.na(country)){
      if(!na_uk){
        if(!na_de){
          return(NA)
        }
        else{
          country <- "UK"
        }
      }  
      else{
        if(!na_de){
          country <- "DE"
        }
        else{
          return(NA)
        }
      }
      #messagef("Fixed country to %s", country)
      if(is.na(country)){
        browser()
      }
    }
    
    #browser()
    if(na_uk && na_de){
      if(is.na(age)){
        return(NA)
      }
      year_group_uk <- year_group_int_from_age(age, "UK")
      year_group_de <- year_group_int_from_age(age, "DE")
      na_uk <- FALSE
      na_de <- FALSE
      if(is.na(year_group_de) || is.na(year_group_uk)){
        browser()
      }
      if((year_group_uk - year_group_de) != 1){
        if(country == "DE"){
          year_group_uk <- NA
          na_uk <- TRUE
        }
        else{
          year_group_de <- NA
          na_de <- TRUE
          
        }
      }
    } 
    if(country == "DE"){
      if(na_de){
        year_group_de <- year_group_uk - 1
      }
      else{
        if(na_uk){
          year_group_uk <- year_group_de + 1
        }
      }
    }
    else if(country == "UK"){
      if(na_uk){
        year_group_uk <- year_group_de + 1
      }
      else{
        if(na_de){
          year_group_de <- year_group_uk - 1
        }
      }
    }
    na_uk <- is.na(year_group_uk)
    na_de <- is.na(year_group_de)
    if(na_uk && na_de){
      return(NA)
    }
    if((year_group_uk - year_group_de) != 1){
      #browser()
      #messagef("Invalid year groups")
      return(NA)
    }
    if(year_group_de == 5 || year_group_de == 13){
      ret <- sprintf("Klasse %02d", year_group_de ) 
    }
    else{
      ret <- sprintf("Klasse %02d/Year %02d", year_group_de, year_group_uk)
    }
    #messagef("-> %s", ret)
    ret
    
  }
)
#' Combined variables are built from sets of numeric variables as their mean with or withou scaling the variables first 
#' For error variables, or variance-type variables, no scaling cam be applied and the aggregate is build using 
#' error propagation rules
#' 
get_combined_var <- function(data, comb_vars, label = "comb_var", scale = T, is_error_var = F){
  l <- 1:nrow(data)
  #browser()
  if(is_error_var && scale){
    stop("z-transform cannot be applied to error vars")
  }
  #data[, comb_vars] <- scale(data[, comb_vars])
  values <- map_dbl(l, function(k){
    if(is_error_var){
      suppressWarnings(sqrt(mean(t(data[k, comb_vars]^2), na.rm = T)))
    }
    else{
      suppressWarnings(mean(t(data[k, comb_vars]), na.rm = T))      
    }
  })
  
  if(scale){
    data[, label] <- scale(values) %>% as.vector()
  }
  else{
    data[, label] <- values
  }
  
  data
}

default_level_vars <- c("level_MIQ" = "MIQ.score",
                        "level_concurrent_music_activity" = "CCM",
                        "level_physical_activity" = "PAC",
                        "level_drama_activity" = "DAC",
                        "level_academic" = "academic.overall",
                        "level_musical_training" = "GMS.musical_training",
                        "level_GMS_general" = "GMS.general",
                        "level_school_difficulties" = "SDQ.difficulties",
                        "level_music_perception" = "music_perception",
                        "level_conscientiousness" = "TPI.conscientiousness",
                        "level_music_home_environment" = "MHE.general_score",
                        "level_cognitive_engagement" = "SEM.cognitive_engagement",
                        "level_behaviorial_engagement" = "SEM.behavioral_engagement",
                        "level_emotional_engagement" = "SEM.emotional_engagement",
                        "level_social_self_concept" = "SCS",
                        "level_academic_self_concept" = "SCA.score",
                        "level_school_difficulties"  = "SDQ.difficulties",
                        "level_prosocial" = "SDQ.prosocial",
                        "level_SES_educational_degree" = "SES.educational_degree")
#' Levelled variables are grouped versiond of some central measurements
#' Grouing is done with the ntile which tries to achieve equal size groups
#' Except leve_drama_activity, all levelled vars have three levesl (Low, Mid, High)
#' 
add_levelled_vars <- function(data, vars = default_level_vars, num_levels = 3, labels = c("Low", "Mid", "High"), aggregate = F){
  common_vars <- intersect(names(data), names(vars))
  data <- data %>% select(-all_of(common_vars))
  if(is.null(names(vars))){
    names(vars) <- vars
  }
  if(aggregate){
    tmp_data <- data %>% group_by(p_id) %>% summarise_at(as.vector(vars), mean, na.rm = T) %>% ungroup()
  }
  else{
    tmp_data <- data
    
  }
  map(1:length(vars), function(i) {
    messagef("Adding %s", names(vars)[i])
    
    if(vars[i] == "level_drama_activity"){
      tmp_data <<- tmp_data %>% mutate(!!names(vars)[i] := factor(ntile(tmp_data %>% pull(vars[i]), 2), labels = labels))
    }
    else{
      tmp_data <<- tmp_data %>% mutate(!!(names(vars)[i]) := factor(ntile(tmp_data %>% pull(vars[i]), num_levels), labels = labels))
 
    }

  })
  if(aggregate){

    tmp_data <- data %>% left_join(tmp_data %>% select(p_id, names(vars)), by = "p_id")
    
  }
  
  as_tibble(tmp_data)
}