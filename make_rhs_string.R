# Writen by Pedro Neves on 06/02/17, under the GPL-2 license.
# Code adapted from package DAISIE (Etienne, Valente, Phillimore & Haegeman), 
# requires package odeintr (Keitt)

library(odeintr)
library(DAISIE)
library(beepr)

#### Testing conditions
lac = pars1[1]
mu = pars1[2]
K = pars1[3]
gam = pars1[4]
laa = pars1[5]
kk = 0
ddep = 0
x <- c(1:11)
lx = (length(x) - 1)/2

nn = -2:(lx+2*kk+1)
lnn = length(nn)          
nn = pmax(rep(0,lnn),nn)

# ddep = 0
laavec = laa * rep(1,lnn)
lacvec = lac * rep(1,lnn)
muvec = mu * rep(1,lnn)
gamvec = gam * rep(1,lnn)


xx1 = c(0,0,x[1:lx],0)
xx2 = c(0,0,x[(lx + 1):(2 * lx)],0)
xx3 = x[2 * lx + 1]

nil2lx = 3:(lx + 2)

il1 = nil2lx+kk-1
il2 = nil2lx+kk+1
il3 = nil2lx+kk
il4 = nil2lx+kk-2

in1 = nil2lx+2*kk-1
in2 = nil2lx+1
in3 = nil2lx+kk

ix1 = nil2lx-1
ix2 = nil2lx+1
ix3 = nil2lx
ix4 = nil2lx-2

# List of pars and indices

list_pars <- list(laavec = laavec, lacvec = lacvec, 
                  muvec = muvec, gamvec = gamvec,
                  nn = nn, xx1 = xx1, xx2 = xx2, xx3 = xx3)
list_indices <- list(il1 = il1, il2 = il2, il3 = il3, il4 = il4,
                     in1 = in1, in2 = in2, in3 = in3,
                     ix1 = ix1, ix2 = ix2, ix3 = ix3, ix4 = ix4)

##### Function to get str of all rhs and pars ####
make_rhs_1 <- function(list_pars, list_indices)
{
  
  # Aux objects
  list_dx <- list()
  pars_list <- list()
  par_name_list <- list()
  init_state_list <- list()
  init_state_list_names <-list()
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
      assign(temp_xx2_ix4, list_pars$xx2[ix4][i])
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
    assign(temp_lacvec_il4_plusone, list_pars$lacvec[il4 + 1][i])
    
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
    
    prod5 <- paste0(" (-1.0 * (", paste(temp_laavec_il3_plusone, temp_xx2_ix3, sep = " * "), "))")
    
    # dx2 rhs of equation
    complete_rhs <- paste(prod1, prod2, prod3, prod4, prod5, sep = " + ")
    list_dx[[lx + dx_list_counter_increase + 1]] <- complete_rhs

    
    #### dx3 ####
    if(i == length(list_pars$laavec[list_indices$il1])){
      
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


make_sys <- function(rhs)
{
  ode_system <- list()
  for(i in 1:length(rhs$rhs)){
    ode_system[[i]] <- (paste0("dxdt[", i - 1, "] = ", rhs$rhs[[i]], "; "))
  }
  
  return(paste(ode_system, sep = "", collapse = ""))
}


#### Run code ####


sys_5_system <- make_rhs_1(list_pars, list_indices)
eqs_5_system <- make_sys(sys_5_system)

sys_101_system <- make_rhs_1(list_pars, list_indices)
eqs_101_system <- make_sys(sys_101_system)

pars_5_system <- sys_5_system$pars
pars_5_system <- pars_5_system[unique(names(pars_5_system))]

pars_101_system <- sys_101_system$pars
pars_101_system <- pars_101_system[unique(names(pars_101_system))]


y_5_system <- compile_sys(name = "y_5_system", eqs_5_system, pars_5_system, 
                          sys_dim = 5, atol = abstol, rtol = reltol) 
beep(sound = 2)


y_101_system <- compile_sys(name = "y_101_system", eqs_101_system, pars_101_system,
                            sys_dim = 101, atol = abstol, rtol = reltol) 
beep(sound = 2)

compile_sys(name = "y", sys$rhs, pars[unique(names(pars))])

x <- c(1,0,0,0,0)
probs_5_system <- x
result_5_system <- y_5_system(x, 4, .10)
deSolve_5_system <- ode(probs_5_system,brts[1:2],DAISIE_loglik_rhs,c(pars1,k1,ddep),
    rtol = reltol,atol = abstol,method = methode)
write.csv(result_5_system, "odeintr_5_system.csv")
write.csv(deSolve_5_system, "deSolve_5_system.csv")

x <- rep(0, 101)
x[1] <- 1
probs_101_system <- x
result_101_system <- y_101_system(x, 4, .10)
deSolve_101_system <- ode(probs_101_system,brts[1:2],DAISIE_loglik_rhs,c(pars1,k1,ddep),
                          rtol = reltol,atol = abstol,method = methode)
write.csv(result_101_system, "odeintr_101_system.csv")
write.csv(deSolve_101_system, "desolve_101_system.csv")



write(y_5_system, "5_component_code.cpp")
write(y_101_system, "101_component_code.cpp")

write(eqs_5_system, "rhs_5.txt")
write(eqs_101_system, "rhs_101.txt")


