# Run code

source("make_rhs_string.R")

ddep_vec <- c(1, 11, 2)
time_to_run <- c()
for (i in 1){
  t_start <- Sys.time()
  test_integrators(300, ddep = ddep_vec[i], K = 40)
  t_end <- Sys.time()
  time_to_run[i] <- t_end - t_start
}

