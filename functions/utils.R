printf <- function(...) print(sprintf(...))
messagef <- function(...) message(sprintf(...))

pull_unique <- function(data, varname){
  data %>% pull(varname) %>% uniqe()  
}