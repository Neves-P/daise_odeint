# Written by Pedro Neves on 06/02/17, under the GPL-2 license.
# Code adapted from package DAISIE (Etienne, Valente, Phillimore & Haegeman), 
# requires package odeintr (Keitt)

prepare_odeintr <- function(probs, pars){
  # Parameter testing function #
  lac <- pars[1]
  mu <- pars[2]
  K <- pars[3]
  gam <- pars[4]
  laa <- pars[5]
  kk <- pars[6]
  ddep <- pars[7]
  x <- probs
  lx <- (length(x) - 1) / 2

  nn <- -2:(lx + 2 * kk + 1)
  lnn <- length(nn)
  nn <- pmax(rep(0, lnn), nn)

  if (ddep == 0)
  {
    laavec = laa * rep(1, lnn)
    lacvec = lac * rep(1, lnn)
    muvec = mu * rep(1, lnn)
    gamvec = gam * rep(1, lnn)
  } else {
    if (ddep == 1)
    {
      laavec = laa * rep(1, lnn)
      lacvec = pmax(rep(0, lnn),lac * (1 - nn / K))
      muvec = mu * rep(1, lnn)
      gamvec = gam * rep(1, lnn)
    } else {
      if (ddep == 2)
      {
        laavec = laa * rep(1, lnn)
        lacvec = pmax(rep(0, lnn),lac * exp(-nn / K))
        muvec = mu * rep(1, lnn)
        gamvec = gam * rep(1, lnn)
      } else {
        if (ddep == 11)
        {
          laavec = laa * rep(1 , lnn)
          lacvec = pmax(rep(0 , lnn), lac * (1 - nn / K))
          muvec = mu * rep(1, lnn)
          gamvec = pmax(rep(0, lnn), gam * (1 - nn / K))
        } else {
          if (ddep == 21)
          {
            laavec = laa * rep(1, lnn)
            lacvec = pmax(rep(0, lnn), lac * exp(-nn / K))
            muvec = mu * rep(1, lnn)
            gamvec = pmax(rep(0, lnn), gam * exp(-nn / K))
          } else {
            if (ddep == 3)
            {
              laavec = laa * rep(1, lnn)
              lacvec = lac * rep(1, lnn)
              muvec = mu * (1 + nn / K)
              gamvec = gam * rep(1, lnn)
            }
          }
        }
      }
    }
  }



  xx1 <- c(0, 0, x[1:lx],0)
  xx2 <- c(0, 0, x[(lx + 1):(2 * lx)], 0)
  xx3 <- x[2 * lx + 1]

  nil2lx <- 3:(lx + 2)

  il1 <- nil2lx + kk - 1
  il2 <- nil2lx + kk + 1
  il3 <- nil2lx + kk
  il4 <- nil2lx + kk - 2

  in1 <- nil2lx + 2 * kk - 1
  in2 <- nil2lx + 1
  in3 <- nil2lx + kk

  ix1 <- nil2lx - 1
  ix2 <- nil2lx + 1
  ix3 <- nil2lx
  ix4 <- nil2lx - 2

  # List of pars and indices

  list_pars <- list(laavec = laavec, lacvec = lacvec,
                    muvec = muvec, gamvec = gamvec,
                    nn = nn, xx1 = xx1, xx2 = xx2, xx3 = xx3, lx = lx)
  list_indices <- list(il1 = il1, il2 = il2, il3 = il3, il4 = il4,
                       in1 = in1, in2 = in2, in3 = in3,
                       ix1 = ix1, ix2 = ix2, ix3 = ix3, ix4 = ix4)

  return(list(list_pars, list_indices))
}

compile_DAISIE <- function(probs, pars, beep){
  # Compiles system from DAISIE in odeintr #

  list_pars_indices <- prepare_odeintr(probs, pars)
  sys <- make_rhs_1(list_pars_indices[[1]], list_pars_indices[[2]])
  eqs <- make_sys(sys)

  pars <- sys$pars[unique(names(sys$pars))]

  compile_sys(name = "y_odeintr", eqs, pars, 
              sys_dim = length(sys$rhs), atol = 1e-10, rtol = 1e-10, 
              rebuild = FALSE, cleanupCacheDir = TRUE) 
  if (beep == TRUE){
  beep(2)
  }
}

integrate_daisie <- function(probs, pars, t = 4, timestep = 0.5){
  # Integrates DAISIE sytem with odeintr and deSolve #

  brts <- c(-t, 0)
  result_deSolve <- ode(probs,
                        brts[1:2], DAISIE_loglik_rhs, pars,
                        rtol = 1e-10, atol = 1e-10, method = "lsodes")
   try(write.csv(result_deSolve, "deSolve.csv"), outFile = "Could not write
       to deSolve.csv.\n")
  
  result_odeintr <- y_odeintr(probs, t, timestep) ### TEST: ODEINTR CRASHES WHEN X IS INVALID
   try(write.csv(result_odeintr, "odeintr.csv"), outFile = "Could not write
       to odeintr.csv.\n")
  
  # Unload DLLs
  # unloadNamespace(y_odeintr)
  # dyn.unload()
  # dir(paste0(tempdir(), "\\sourceCpp-x86_64-w64-mingw32-0.12.15")tempdir(), pattern = "sourceCpp")
  return(list(deSolve = result_deSolve, odeintr = result_odeintr))
}

##### Function to get str of all rhs and pars ####
make_rhs_1 <- function(list_pars, list_indices){

  # Aux objects
  lx <- list_pars$lx
  list_dx <- list()
  pars_list <- list()
  par_name_list <- list()
  dx_list_counter <- 1
  dx_list_counter_increase <- 0
  x_counter <- 0
  x_counter_2 <- lx


  for (i in 1:length(list_pars$laavec[list_indices$il1])){
    
    # Generate X first
    
    ## xx1 = Qk
    
    if(i != 1){ # n-1 falls out of boundary
      temp_xx1_ix1 <- paste0("x[", x_counter - 1, "]")
      assign(temp_xx1_ix1, list_pars$xx1[list_indices$ix1][i])
    }else{
      temp_xx1_ix1 <- "0.0"
    }
    
    if(i != length(list_pars$laavec[list_indices$il1])){ # n+1 falls out of boundary
      temp_xx1_ix2 <- paste0("x[", x_counter + 1, "]")
      assign(temp_xx1_ix2, list_pars$xx1[list_indices$ix2][i])
    }else{
      temp_xx1_ix2 <- "0.0"
    }
    
    ### Qk,n (n term, no -1 or +1 needed)
    temp_xx1_ix3 <- paste0("x[", x_counter, "]")
    assign(temp_xx1_ix3, list_pars$xx1[list_indices$ix3][i])
    
    ## xx2 = QMk
    
    if(i != 1){ # n-1 falls out of boundary
      temp_xx2_ix1 <- paste0("x[", x_counter_2 - 1, "]")
      assign(temp_xx2_ix1, list_pars$xx2[list_indices$ix1][i])
    }else{
      temp_xx2_ix1 <- "0.0"
    }
    
    if(i != length(list_pars$laavec[list_indices$il1])){ # n+1 falls out of boundary
      temp_xx2_ix2 <- paste0("x[", x_counter_2 + 1, "]")
      assign(temp_xx2_ix2, list_pars$xx2[list_indices$ix2][i])
    }else{
      temp_xx2_ix2 <- "0.0"
    }
    
    ### QMk,n (n term, no -1 or +1 needed)
      temp_xx2_ix3 <- paste0("x[", x_counter_2, "]")
      assign(temp_xx2_ix3, list_pars$xx2[list_indices$ix3][i])
    
    if(i > 2 ){ # n - 2 falls out of boundary
      temp_xx2_ix4 <- paste0("x[", x_counter_2 - 2, "]")
      assign(temp_xx2_ix4, list_pars$xx2[list_indices$ix4][i])
    }else{
      temp_xx2_ix4 <- "0.0"
    }


    #### dx1 ####

    # First product
    temp_laavec_il1_plusone <- paste("laavec_il1_plusone", i, sep = "_")
    assign(temp_laavec_il1_plusone, list_pars$laavec[list_indices$il1 + 1][i])


    prod1 <- paste(temp_laavec_il1_plusone, temp_xx2_ix1, sep = " * ")


    # Second product
    temp_lacvec_il4_plusone <- paste("lacvec_il4_plusone", i, sep = "_")
    assign(temp_lacvec_il4_plusone, list_pars$lacvec[list_indices$il4 + 1][i])

    prod2 <- paste(temp_lacvec_il4_plusone, temp_xx2_ix4, sep = " * ")


    # Third product
    temp_muvec_il2_plusone <- paste("muvec_il2_plusone", i, sep = "_")
    assign(temp_muvec_il2_plusone, list_pars$muvec[list_indices$il2 + 1][i])

    prod3 <- paste(temp_muvec_il2_plusone, temp_xx2_ix3, sep = " * ")

    # Fourth product
    temp_lacvec_il1 <- paste("lacvec_il1", i, sep = "_")
    assign(temp_lacvec_il1, list_pars$lacvec[list_indices$il1][i])

    temp_nn_in1 <- paste("nn_in1", i, sep = "_")
    assign(temp_nn_in1, list_pars$nn[list_indices$in1][i])

    prod4 <- paste(temp_lacvec_il1, temp_nn_in1, temp_xx1_ix1, sep = " * ")


    # Fifth product
    temp_muvec_il2 <- paste("muvec_il2", i, sep = "_")
    assign(temp_muvec_il2, list_pars$muvec[list_indices$il2][i])

    temp_nn_in2 <- paste("nn_in2", i, sep = "_")
    assign(temp_nn_in2, list_pars$nn[list_indices$in2][i])

    prod5 <- paste(temp_muvec_il2, temp_nn_in2, temp_xx1_ix2, sep = " * ")

    # Negative term
    temp_muvec_il3 <- paste("muvec_il3", i, sep= "_")
    assign(temp_muvec_il3, list_pars$muvec[list_indices$il3][i])

    temp_lacvec_il3 <- paste("lacvec_il3", i, sep = "_")
    assign(temp_lacvec_il3, list_pars$lacvec[list_indices$il3][i])

    neg_term1 <- paste("(-1.0) * (", paste(temp_muvec_il3, temp_lacvec_il3, 
                                           sep = " + "), ")", sep = "")


    # Sixth product
    temp_nn_in3 <- paste("nn_in3", i, sep = "_")
    assign(temp_nn_in3, list_pars$nn[list_indices$in3][i])

    prod6 <- paste(neg_term1, temp_nn_in3, temp_xx1_ix3, sep = " * ")

    # Seventh product
    # This gamvec is assign the regular (non negative) value of gam
    # The string is built with a negative sign for correct input in odeintr
    temp_neggamvec_il3 <- paste("gamvec_il3", i, sep = "_")
    assign(temp_neggamvec_il3, list_pars$gamvec[list_indices$il3][i]) 

    prod7 <- paste(paste0("(-1.0) * ", temp_neggamvec_il3),
                   temp_xx1_ix3, sep = " * ")


    # dx1 rhs of equation
    complete_rhs <- paste(prod1, prod2, prod3, prod4,
                          prod5, prod6, prod7, sep = " + ")
    list_dx[[dx_list_counter + dx_list_counter_increase]] <- complete_rhs


    #### dx2 ####

    # First product
    temp_gamvec_il3 <- paste("gamvec_il3", i, sep = "_")
    assign(temp_gamvec_il3, list_pars$gamvec[list_indices$il3][i])

    prod1 <- paste(temp_gamvec_il3, temp_xx1_ix3, sep = " * ")

    # Second product
    temp_lacvec_il1_plusone <- paste("lacvec_il1_plusone", i, sep = "_")
    assign(temp_lacvec_il1_plusone, list_pars$lacvec[list_indices$il1 + 1][i])

    temp_nn_in1 <- paste("nn_in1", i, sep = "_")
    assign(temp_nn_in1, list_pars$nn[list_indices$in1][i])

    prod2 <- paste(temp_lacvec_il1_plusone, temp_nn_in1,
                   temp_xx2_ix1, sep = " * ")

    # Third product
    temp_muvec_il2_plusone <- paste("muvec_il2_plusone", i, sep = "_")
    assign(temp_muvec_il2_plusone, list_pars$muvec[list_indices$il2 + 1][i])

    temp_nn2_in2 <- paste("nn_in2", i, sep = "_")
    assign(temp_nn2_in2, list_pars$nn[list_indices$in2][i])
    
    prod3 <- paste(temp_muvec_il2_plusone,
                   temp_nn_in2, temp_xx2_ix2, sep = " * ")
    
    # Negative term
    temp_muvec_il3_plusone <- paste("muvec_il3_plusone", i, sep= "_")
    assign(temp_muvec_il3_plusone, list_pars$muvec[list_indices$il3 + 1][i])
    
    temp_lacvec_il3_plusone <- paste("lacvec_il3_plusone", i, sep = "_")
    assign(temp_lacvec_il3_plusone, list_pars$lacvec[list_indices$il3 + 1][i])
    
    neg_term1 <- paste("-(", paste(temp_muvec_il3_plusone, 
                                   temp_lacvec_il3_plusone, 
                                   sep = " + "), ")", sep = "")
    
    # Fourth product
    temp_nn_in3_plusone <- paste("nn_in3_plusone", i, sep = "_")
    assign(temp_nn_in3_plusone, list_pars$nn[list_indices$in3 + 1][i])
    
    prod4 <- paste(neg_term1, temp_nn_in3_plusone, temp_xx2_ix3, sep = " * ")
    
    # Fifth product
    temp_laavec_il3_plusone <- paste("laavec_il3_plusone", i, sep = "_")
    assign(temp_laavec_il3_plusone, list_pars$laavec[list_indices$il3 + 1][i])
    
    prod5 <- paste0(" (-1.0 * (", paste(temp_laavec_il3_plusone,
                                        temp_xx2_ix3, sep = " * "), "))")
    
    # dx2 rhs of equation
    complete_rhs <- paste(prod1, prod2, prod3, prod4, prod5, sep = " + ")
    list_dx[[lx + dx_list_counter_increase + 1]] <- complete_rhs

    
    #### dx3 ####
    if (i == length(list_pars$laavec[list_indices$il1])){
      
      # Negative term
      temp_laavec_il3_one <- "laavec_il3_one"
      assign(temp_laavec_il3_one, list_pars$laavec[list_indices$il3][1])
      
      temp_lacvec_il3_one <- "lacvec_il3_one"
      assign(temp_lacvec_il3_one, list_pars$lacvec[list_indices$il3][1])
      
      temp_gamvec_il3_one <- "gamvec_il3_one"
      assign(temp_gamvec_il3_one, list_pars$gamvec[list_indices$il3][1])
      
      temp_muvec_il3_one <- "muvec_il3_one"
      assign(temp_muvec_il3_one, list_pars$muvec[list_indices$il3[1]])
      
      neg_term1 <- paste("-(", paste(temp_laavec_il3_one, temp_lacvec_il3_one,
                                     temp_gamvec_il3_one, temp_muvec_il3_one,
                                     sep = " + "), ")", sep = "")
      
      temp_xx3 <- paste0("x[", 2 * lx, "]")
      assign(temp_xx3, list_pars$xx3)
      prod1 <- paste(neg_term1, temp_xx3, sep = " * ")
      
      list_dx[[lx + dx_list_counter_increase + 2]] <- prod1
    }
    
    #### Model parameters per rhs ####
    
    # Updates index of dx_list for next equation loop
    x_counter <- x_counter + 1
    x_counter_2 <- x_counter_2 + 1
    dx_list_counter_increase <- dx_list_counter_increase + 1
    
    # Store initial state and parameters
    local_env_pars <- ls(pattern = "temp")
    
    pars <- local_env_pars[grepl("temp", local_env_pars) &
                             !grepl("xx", local_env_pars)]
    par_name_list[[i]] <- unlist(unname(mget(pars)))
    pars_list[[i]] <- mget(unlist(unname(mget(pars))))
    
  }  
  
  # Return
  
  return(list(rhs = list_dx, pars = unlist(pars_list)))
}

# Appends lhs of equation and ; to system for odeintr to interpret
make_sys <- function(rhs){
  ode_system <- list()
  for (i in 1:length(rhs$rhs)){
    ode_system[[i]] <- (paste0("dxdt[", i - 1, "] = ", rhs$rhs[[i]], "; "))
  }
  
  return(paste(ode_system, sep = "", collapse = ""))
}


#### Run code ####

# Parameter testing


# Generates list of increasingly large probability vectors
# Only keeps vectors of odd size
make_prob_test_list <- function(nruns, start_size = 5){
  if (start_size < 5){
    stop("Start size must be > 5.")
  }
  
  probs_test_list <- list()
  for (i in 1:(nruns * 2)){
    
    # Starts start_size at 5 and increases 1 per loop until i = nruns
    probs_test <- rep(0, start_size - (start_size - (5 + (i - 1))))
    if (length(probs_test) > 201){
      probs_test <- rep(0,201)
    }
    probs_test[1] <- 1
    probs_test_list[[i]] <- probs_test
  }
  index_vector <- c(TRUE,FALSE)
  probs_test_list <- probs_test_list[index_vector]
  return(probs_test_list)
}

# Generates list of parameter vectors based on sampling from a normal 
# distribution
make_pars_test_list <- function(nruns, K = Inf, ddep = 0, kk = 0, seed) {
  set.seed(seed)
  pars_test_list <- list()
  for (i in 1:nruns){
    pars_test_list[[i]] <- c(lac = abs(rnorm(1, 2.5, 1)), 
                             mu = abs(rnorm(1, 2.7, 1)),
                             K = K, gam = abs(rnorm(1, 0.009, 0.05)),
                             laa = abs(rnorm(1, 1.01, 1)), kk = 0, ddep = ddep)
  }
  return(pars_test_list)
}

# Compiles odeintr integrator and integrates system with odeintr and deSolve
run_integrator_test <- function(probs, pars, nruns, beep) {
  pb <- txtProgressBar(min = 0, max = nruns, initial = 1, style = 3)
  
  result_list <- list()
  for (i in 1:nruns){
    
    if (i == 1 && (file.exists(paste0("results_deSolve_ddep", pars[[i]][[7]],
                                      ".csv")) ||
                   file.exists(paste0("results_odeintr_ddep", pars[[i]][[7]],
                                      ".csv")))){
      file.remove(paste0("results_deSolve_ddep", pars[[i]][[7]],
                         ".csv"))
      file.remove(paste0("results_odeintr_ddep", pars[[i]][[7]],
                         ".csv"))
    }
    
    cat(paste0("\nIntegrating system ", i, "...",  "\n"))
    cat(probs[[i]], file="probs.csv", append = FALSE, sep = "\n")
    cat(pars[[i]], file="pars.csv", append = FALSE, sep = "\n")
    compile_DAISIE(probs[[i]], pars[[i]], beep)
    result_list[[i]] <- integrate_daisie(probs = probs[[i]],
                                         pars = pars[[i]],
                                         t = 4, timestep = 0.5)
    setTxtProgressBar(pb, i)
    
    ### Unload DLLs ###
    # Get all DLLs
    loaded_dll <- getLoadedDLLs()
    
    # Get only sourceCpp DLL object
    loaded_dll_cpp <- loaded_dll[grep("sourceCpp", loaded_dll)]
    
    # Path to loaded sourceCpp DLL
    path <- loaded_dll_cpp[1][[1]][[2]]
    
    # Unload DLL
    dyn.unload(path[[1]])
    
    # Write results to file
    write.table(result_list[[i]]$deSolve, 
                file = paste0("results_deSolve_", pars[[i]][[7]], ".csv"),
                append = TRUE, sep = ",")
    write.table(result_list[[i]]$odeintr,
                file = paste0("results_odeintr_", pars[[i]][[7]], ".csv"),
                append = TRUE, sep = ",")
    
    Sys.sleep(0.5)
  }
  return(result_list)
}


# Calculates difference between deSolve and odeintr last equation
calculate_error <- function(results, nruns) {
  error <- c()
  for (i in 1:nruns){
    deSolve_last_col <- results[[i]][[1]][[2,ncol(results[[i]][[1]]) - 1]]
    odeintr_last_col <- results[[i]][[2]][length(results[[i]][[2]]) - 1][nrow(results[[i]][[2]][length(results[[i]][[2]]) - 1]),]
    
    error[i] <- deSolve_last_col - odeintr_last_col
  }
  return(list(error))
}

# Plots error as scatterplot or boxplot
plot_error <- function(error, K, ddep, type = "boxplot") {
  if (type != "boxplot" && type != "scatterplot"){
    cat("\nInvalid plot type. Coercing to boxplot.\n")
    type = "boxplot"
  }
  if (type == "boxplot"){
    boxplot(error, main = "Differences between deSolve and odeintr", 
            ylab = "deSolve last component - odeintr last component")
    legend(x = "bottom", paste0("K = ", K, " --- ddep = ", ddep),
           inset = c(0,-0.2), xpd = TRUE) 
  }else{
    plot(error[[1]], main = "Differences between deSolve and odeintr",
         ylab = "deSolve last component - odeintr last component",
         xlab = "System")
    abline(0,0, col = "red")
    legend(x = "bottom", paste0("K = ", K, " --- ddep = ", ddep),
           inset = c(0,-0.2), xpd = TRUE) 
  }
}

# Compiles and integrates two systems with random parameters and specified size
# Returns results of integration and difference between second to last 
# components of the system.
test_integrators <- function(nruns, start_size = 5, ddep = 0, kk = 0,
                             K = Inf, seed = 42, beep = FALSE, 
                             plotit = TRUE, style = "boxplot") {
  require(DAISIE)
  require(deSolve)
  require(odeintr)
  require(beepr)
  
  # Generate probabilities vector c(1, 0, 0, 0, 0) with odd number of elements
  probs <- make_prob_test_list(nruns, start_size)
  
  # Generate parameters vector from rnorm
  test_pars <- make_pars_test_list(nruns, K, ddep, kk, seed)
  
  # Integrates generated systems using deSolve and odeintr
  results <- run_integrator_test(probs, test_pars, nruns, beep)
  
  # Calculates difference between outputs
  error <- calculate_error(results, nruns)
  
  # Plots differences between outputs
  if (plotit == TRUE){
    plot_error(error, K, ddep, style)
  }
  invisible(list(results, error))
}


#### Benchmarking ####

deSolve_calc <- function(x, pars, brts, timestep = 0.5) {
  return(ode(x,
             brts[1:2], DAISIE_loglik_rhs, pars,
             rtol = 1e-10, atol = 1e-16, method = "lsodes"))
}

odeintr_calc <- function(x, t, timestep){
  y_odeintr(x, t, timestep)
}

