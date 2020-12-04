printf <- function(...) print(sprintf(...))
messagef <- function(...) message(sprintf(...))

pull_unique <- function(data, varname){
  data %>% pull(varname) %>% uniqe()  
}

limiter <- function(x, limits){
  if(is.null(limits) || length(limits) != 2){
    return(x)
  }
  pmax(pmin(x, limits[2]), limits[1])
}

se <- function(x,...){
  sd(x, ...)/sqrt(length(x))
}
