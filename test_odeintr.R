# Writen by Pedro Neves on 13/12/, under the GPL-2 license.
# Code adapted from package DAISIE (Etienne, Valent, Phillimore & Haegeman), 
# requires package odeintr (Keitt)
# 
# 
# RHS1_odeintr --------------------------------------------------------------------


DAISIE_loglik_rhs_odeintr = function(t,x,pars)
{
  # Returns C++ string with rhs of the equation and parameter vector to pass to 
  # odeintr. This code adapted from functions DAISIE_loglik_rhs and 
  # DAISIE_loglik_rhs1, in function DAISIE_loglik_all
  
  # Reads model parameters
  lx = (length(x) - 1)/2 # Number of variables containing the state of the system . x = current probs
  lac = pars[1]          # Lambda^c
  mu = pars[2]           # mu
  K = pars[3]            # K
  gam = pars[4]          # gamma
  laa = pars[5]          # Lambda^a
  kk = pars[6]           # k (nº species in )
  ddep = pars[7]         # Type of diversity dependence (0 = no DD; 1 = linear 
                         # dep in speciation; 2 = exponential dep in 
                         # speciation rate; 11 = linear dep in speciation and 
                         # immigration; 21 = exponential dep in speciation and immigration) 
  
  nn = -2:(lx+2*kk+1)
  lnn = length(nn)         # 
  nn = pmax(rep(0,lnn),nn) # Number of species in clade
  
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
          }}}}}
  
  #x = x * (x > 0)
  
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

  # Creates vector with parameters needed for ODE integration. Note: order of 
  # these parameters matters, as they are assigned their name to be called 
  # below. TODO(Neves-P): should make this dynamic and not "magic number"ish
  model.pars <- c(
    laavec[il1 + 1],
    laavec[il3[1]],
    laavec[il3 + 1],
    laavec[il3],
    lacvec[il4 + 1],
    lacvec[il1],
    lacvec[il3],
    lacvec[il3[1]],
    lacvec[il1 + 1],
    lacvec[il3 + 1],
    xx1[ix1],
    xx1[ix2],
    xx1[ix3],
    xx2[ix1],
    xx2[ix2],
    xx2[ix3],
    xx2[ix4],
    xx3,
    nn[in1],
    nn[in2],
    nn[in3],
    muvec[il2 + 1],
    muvec[il2],
    muvec[il3],
    muvec[il3 + 1],
    muvec[il3[1]],
    gamvec[il3],
    gamvec[il3[1]])
  
  # Assigns C++ friendly parameter names to model parameters
  model.pars.names <- c(
    "laa_il1_add_one",
    "laa_l3_1",
    "laa_il3_add_one",
    "laa_il3",
    "lac_il4_add_one",
    "lac_il1",
    "lac_il3",
    "lac_il3_1",
    "lac_il1_add_one",
    "lac_il3_add_one",
    "xx1_ix1",
    "xx1_ix2",
    "xx1_ix3",
    "xx2_ix1",
    "xx2_ix2",
    "xx2_ix3",
    "xx2_ix4",
    "xx3",
    "nn_in1",
    "nn_in2",
    "nn_in3",
    "mu_il2_add_one",
    "mu_il2",
    "mu_il3",
    "mu_il3_add_one",
    "mu_il3_1",
    "gam_il3",
    "gam_il3_1")

  names(model.pars) <- modelpars.names
  

  

# TODO(Neves-P): Do a paste function here that prints the sum of all the elements
#  of each lac[il1]etc. vector. May not be necessary

# Stores C++ code in string to be passed to compile_sys function. Arguments are
# passed by their C++ friendly name (names(modelpars)).
# 
rhs1.cpp.system <- 'dxdt[0] = laa_il1_add_one * xx2_ix1 + lac_il4_add_one * xx2_ix4 + mu_il2_add_one * xx2_ix3 + lac_il1 * nn_in1 * xx1_ix1 + mu_il2 * nn_in2 * xx1_ix2 - (mu_il3 + lac_il3) * nn_in3 * xx1_ix3 + gam_il3 * xx1_ix3; dxdt[1] = gam_il3 * xx1_ix3 + lac_il1_add_one * nn_in1 * xx2_ix1 + mu_il2_add_one * nn_in2 * xx2_ix2 - (mu_il3_add_one + laa_il3_add_one) * nn_in3 * xx2_ix3 - gam_il3 * xx1_ix3; dxdt[2] = -(laa_il3 + lac_il3 + gam_il3_1 + mu_il3_1) * xx3;'

# TODO(Neves-P): Build rhs2 for k=1

return(list(rhs1.cpp.system, model.pars))

# Prints the C++ code built by odeintr. Move to another function?
# the_code <- compile_sys("rhs1" , rhs1.cpp.system, model.pars, TRUE)
}

