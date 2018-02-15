# Written by Pedro Neves on 06/02/17, under the GPL-2 license.
# Code adapted from package DAISIE (Etienne, Valente, Phillimore & Haegeman), 
# requires package odeintr (Keitt)

library(odeintr)
library(DAISIE)
library(beepr)
library(deSolve)

prepare_odeintr <- function(pars, x){
  # Parameter testing function #

  lac <- pars[1]
  mu <- pars[2]
  K <- pars[3]
  gam <- pars[4]
  laa <- pars[5]
  kk <- pars[6]
  ddep <- pars[7]
  x <- x
  lx <- (length(x) - 1) / 2
  
  nn <- -2:(lx + 2 * kk + 1)
  lnn <- length(nn)
  nn <- pmax(rep(0, lnn), nn)

  if(ddep == 0)
  {
    laavec = laa * rep(1,lnn)
    lacvec = lac * rep(1,lnn)
    muvec = mu * rep(1,lnn)
    gamvec = gam * rep(1,lnn)
  } else {
    if(ddep == 1)
    {
      laavec = laa * rep(1,lnn)
      lacvec = pmax(rep(0,lnn),lac * (1 - nn/K))
      muvec = mu * rep(1,lnn)
      gamvec = gam * rep(1,lnn)
    } else {
      if(ddep == 2)
      {
        laavec = laa * rep(1,lnn)
        lacvec = pmax(rep(0,lnn),lac * exp(-nn/K))
        muvec = mu * rep(1,lnn)
        gamvec = gam * rep(1,lnn)
      } else {
        if(ddep == 11)
        {
          laavec = laa * rep(1,lnn)
          lacvec = pmax(rep(0,lnn),lac * (1 - nn/K))
          muvec = mu * rep(1,lnn)
          gamvec = pmax(rep(0,lnn),gam * (1 - nn/K))
        } else {
          if(ddep == 21)
          {
            laavec = laa * rep(1,lnn)
            lacvec = pmax(rep(0,lnn),lac * exp(-nn/K))
            muvec = mu * rep(1,lnn)
            gamvec = pmax(rep(0,lnn),gam * exp(-nn/K))
          } else {
            if(ddep == 3)
            {
              laavec = laa * rep(1,lnn)
              lacvec = lac * rep(1,lnn)
              muvec = mu * (1 + nn/K)
              gamvec = gam * rep(1,lnn)
            }
          }
        }
      }
    }
  }



  xx1 <- c(0,0,x[1:lx],0)
  xx2 <- c(0,0,x[(lx + 1):(2 * lx)],0)
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

compile_DAISIE <- function(pars, x){
  # Compiles system from DAISIE in odeintr #
  
  x <- x
  list_pars_indices <- prepare_odeintr(pars, x)
  sys <- make_rhs_1(list_pars_indices[[1]], list_pars_indices[[2]])
  eqs <- make_sys(sys)
  
  pars <- sys$pars[unique(names(sys$pars))]
  
  compile_sys(name = "y_odeintr", eqs, pars, 
              sys_dim = length(sys$rhs), atol = 1e-10, rtol = 1e-10) 
  beep(2)
}

integrate_DAISIE <- function(x, pars, t = 4, timestep = 0.5){
  # Integrates DAISIE sytem with odeintr and deSolve #
  brts <- c(-t, 0)
  result_deSolve <- ode(x,
                        brts[1:2], DAISIE_loglik_rhs, pars,
                        rtol = 1e-10,atol = 1e-10,method = "lsodes")
  
  result_odeintr <- y_odeintr(x, t, timestep)
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
  
  
  for(i in 1:length(list_pars$laavec[list_indices$il1])){
    
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
    
    
    prod1 <- paste(temp_laavec_il1_plusone, temp_xx2_ix1 , sep = " * ")
    
    
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
    
    prod7 <- paste(paste0("(-1.0) * " , temp_neggamvec_il3), temp_xx1_ix3, sep = " * ")
    
    
    # dx1 rhs of equation
    complete_rhs <- paste(prod1, prod2, prod3, prod4, prod5, prod6, prod7, sep = " + ")
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
    
    prod2 <- paste(temp_lacvec_il1_plusone, temp_nn_in1, temp_xx2_ix1, sep = " * ")
    
    # Third product
    temp_muvec_il2_plusone <- paste("muvec_il2_plusone", i, sep = "_")
    assign(temp_muvec_il2_plusone, list_pars$muvec[list_indices$il2 + 1][i])
    
    temp_nn2_in2 <- paste("nn_in2", i, sep = "_")
    assign(temp_nn2_in2, list_pars$nn[list_indices$in2][i])
    
    prod3 <- paste(temp_muvec_il2_plusone, temp_nn_in2, temp_xx2_ix2, sep = " * ")
    
    # Negative term
    temp_muvec_il3_plusone <- paste("muvec_il3_plusone", i, sep= "_")
    assign(temp_muvec_il3_plusone, list_pars$muvec[list_indices$il3 + 1][i])
    
    temp_lacvec_il3_plusone <- paste("lacvec_il3_plusone", i, sep = "_")
    assign(temp_lacvec_il3_plusone, list_pars$lacvec[list_indices$il3 + 1][i])
    
    neg_term1 <- paste("-(", paste(temp_muvec_il3_plusone, temp_lacvec_il3_plusone, 
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


make_sys <- function(rhs){
  ode_system <- list()
  for (i in 1:length(rhs$rhs)){
    ode_system[[i]] <- (paste0("dxdt[", i - 1, "] = ", rhs$rhs[[i]], "; "))
  }
  
  return(paste(ode_system, sep = "", collapse = ""))
}


#### Run code ####

# Parameter testing

probs_test_list <- list()
for (i in 5:506){
  probs_test_list[[i - 4]] <- rep(0, i)
  probs_test_list[[i - 4]][1] <- 1
}

odd_numbers <-  seq(1, 506, by = 2)
odd_probs <- probs_test_list[odd_numbers]

params_test_list <- list()
for (i in 1:10){
  params_test_list[[i]] <- c(lac = abs(rnorm(1, 2.5, 1)), 
                             mu = abs(rnorm(1, 2.7, 1)),
                             K = Inf, gam = abs(rnorm(1, 0.009, 0.05)),
                             laa = abs(rnorm(1, 1.01, 1)), kk = 0, ddep = 0)
}

result_list <- list()
for (i in 1:10){
  cat(paste0("Integrating function ", i, "...",  "\n"))
  compile_DAISIE(params_test_list[[i]], x = odd_probs[[i]])
  result_list[[i]] <- integrate_DAISIE(x = odd_probs[[i]],
                                       pars = params_test_list[[i]],
                                       t = 4, timestep = 1)
  Sys.sleep(0.5)
}

# Calculates difference between deSolve and odeintr last equation
calculate_error <- function(result_list) {
    diff <- c()
    for (i in 1:10){
    diff <- append(diff, result_list[[i]]$deSolve[length(result_list[[i]]$deSolve) - 2] -
                     result_list[[i]]$odeintr[5,length(result_list[[i]]$odeintr) - 1])
  }
  return(diff)
}

deSolve_calc <- function(x, pars, brts, timestep = 0.5) {
  return(ode(x,
             brts[1:2], DAISIE_loglik_rhs, pars,
             rtol = 1e-10,atol = 1e-16,method = "lsodes"))
}

odeintr_calc <- function(x, t, timestep){
  y_odeintr(x, t, timestep)
}

