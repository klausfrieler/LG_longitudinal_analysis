library(tidyverse)
library(aws.s3)

#source("functions/style.R")

#' Checks whether AWS credentials are already set; if they are not,
#' then sets the credentials.
aws.set_credentials_if_unset <- function(file = "aws-credentials.txt") {
  if (!aws.are_credentials_set()) aws.set_credentials(file) else message("Using pre-loaded AWS credentials.")
}

aws.are_credentials_set <- function() {
  cred <- Sys.getenv(c("AWS_ACCESS_KEY_ID", "AWS_SECRET_ACCESS_KEY"))
  !any(cred == "")
}

aws.clear_credentials <- function() {
  Sys.setenv(
    "AWS_ACCESS_KEY_ID" = "",
    "AWS_SECRET_ACCESS_KEY" = ""
  )
}

aws.set_credentials <- function(file) {
  cred <- if (file.exists(file)) {
    message("Reading AWS credentials from file.")
    aws.get_cred_from_file(file)
  } else {
    message("No credentials found at ", file, ".")
    message("You can avoid entering AWS credentials manually by storing them ",
            "in a text file at this file location.")
    message("The file should comprise two lines of text, plus a final line break.")
    message("The first line should be your AWS access key ID.")
    message("The second line should be your AWS secret access key.")
    aws.get_cred_from_user()
  }
  assertthat::assert_that(
    all(names(cred) == c("id", "secret"))
  )
  message("Setting AWS credentials.")
  Sys.setenv("AWS_ACCESS_KEY_ID" = cred$id,
             "AWS_SECRET_ACCESS_KEY" = cred$secret)
}

aws.get_cred_from_file <- function(file) {
  tryCatch({
    cred <- readLines(file) %>% 
      (function(x) list(id = x[1], secret = x[2]))
    if (is.na(cred$secret)) stop()
  },
  error = function(e) {
    message("Failed to read credentials from ", file, ".")
    message("The file should comprise two lines of text, plus a final line break.")
    message("The first line should be your AWS access key ID.")
    message("The second line should be your AWS secret access key.")
    message("Try fixing ", file, " and rerunning the script, or alternatively delete the file",
            " and enter the AWS credentials manually from the console when prompted.")
    stop("Failed to read AWS credentials from file.")
  })
  cred
}

aws.get_cred_from_user <- function() {
  list(
    id = readline(prompt = "Please enter your AWS access key ID: "),
    secret = readline(prompt = "Please enter your AWS secret access key: ")
  )
}
#if(check_style) style.check_all_functions_with_prefix("aws.")