

# RHS1 --------------------------------------------------------------------

DAISIE_loglik_rhs = function(t,x,pars)
{
	lx = (length(x) - 1)/2
	lac = pars[1]
	mu = pars[2]
	K = pars[3]
	gam = pars[4]
	laa = pars[5]
	kk = pars[6]
	ddep = pars[7]
	
	nn = -2:(lx+2*kk+1)
	lnn = length(nn)
	nn = pmax(rep(0,lnn),nn)
	
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
	
	dx1 = laavec[il1 + 1] * xx2[ix1] + lacvec[il4 + 1] * xx2[ix4] + muvec[il2 + 1] * xx2[ix3] +
		lacvec[il1] * nn[in1] * xx1[ix1] + muvec[il2] * nn[in2] * xx1[ix2] +
		-(muvec[il3] + lacvec[il3]) * nn[in3] * xx1[ix3] +
		-gamvec[il3] * xx1[ix3]
	dx1[1] = dx1[1] + laavec[il3[1]] * xx3 * (kk == 1)
	dx1[2] = dx1[2] + 2 * lacvec[il3[1]] * xx3 * (kk == 1)
	
	dx2 = gamvec[il3] * xx1[ix3] +
		lacvec[il1 + 1] * nn[in1] * xx2[ix1] + muvec[il2 + 1] * nn[in2] * xx2[ix2] +
		-(muvec[il3 + 1] + lacvec[il3 + 1]) * nn[in3 + 1] * xx2[ix3] +
		-laavec[il3 + 1] * xx2[ix3]
	
	dx3 = -(laavec[il3[1]] + lacvec[il3[1]] + gamvec[il3[1]] + muvec[il3[1]]) * xx3
	
	return(list(c(dx1,dx2,dx3)))
}

DAISIE_loglik_rhs2 = function(t,x,pars)
{
	lx = length(x)/3
	lac = pars[1]
	mu = pars[2]
	K = pars[3]
	gam = pars[4]
	laa = pars[5]
	kk = pars[6]
	ddep = pars[7]
	
	nn = -2:(lx+2*kk+1)
	lnn = length(nn)
	nn = pmax(rep(0,lnn),nn)
	
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
	xx3 = c(0,0,x[(2 * lx + 1):(3 * lx)],0)
	
	nil2lx = 3:(lx + 2)
	
	il1 = nil2lx+kk-1
	il2 = nil2lx+kk+1
	il3 = nil2lx+kk
	il4 = nil2lx+kk-2
	
	in1 = nil2lx+2*kk-1
	in2 = nil2lx+1
	in3 = nil2lx+kk
	in4 = nil2lx-1
	
	ix1 = nil2lx-1
	ix2 = nil2lx+1
	ix3 = nil2lx
	ix4 = nil2lx-2
	
	# inflow:
	# anagenesis in colonist when k = 1: Q_M,n -> Q^1_n; n+k species present
	# cladogenesis in colonist when k = 1: Q_M,n-1 -> Q^1_n; n+k-1 species present; rate twice
	# anagenesis of reimmigrant: Q^M,k_n-1 -> Q^k,n; n+k-1+1 species present
	# cladogenesis of reimmigrant: Q^M,k_n-2 -> Q^k,n; n+k-2+1 species present; rate once
	# extinction of reimmigrant: Q^M,k_n -> Q^k,n; n+k+1 species present
	# cladogenesis in one of the n+k-1 species: Q^k_n-1 -> Q^k_n; n+k-1 species present; rate twice for k species
	# extinction in one of the n+1 species: Q^k_n+1 -> Q^k_n; n+k+1 species present
	# outflow:
	# all events with n+k species present
	dx1 = (laavec[il3] * xx3[ix3] + 2 * lacvec[il1] * xx3[ix1]) * (kk == 1) + 
		laavec[il1 + 1] * xx2[ix1] + lacvec[il4 + 1] * xx2[ix4] + muvec[il2 + 1] * xx2[ix3] +
		lacvec[il1] * nn[in1] * xx1[ix1] + muvec[il2] * nn[in2] * xx1[ix2] +
		-(muvec[il3] + lacvec[il3]) * nn[in3] * xx1[ix3] - gamvec[il3] * xx1[ix3]
	
	# inflow:
	# immigration when there are n+k species: Q^k,n -> Q^M,k_n; n+k species present
	# cladogenesis in n+k-1 species: Q^M,k_n-1 -> Q^M,k_n; n+k-1+1 species present; rate twice for k species
	# extinction in n+1 species: Q^M,k_n+1 -> Q^M,k_n; n+k+1+1 species present
	# outflow:
	# all events with n+k+1 species present
	dx2 = gamvec[il3] * xx1[ix3] +
		lacvec[il1 + 1] * nn[in1] * xx2[ix1] + muvec[il2 + 1] * nn[in2] * xx2[ix2] +
		-(muvec[il3 + 1] + lacvec[il3 + 1]) * nn[in3 + 1] * xx2[ix3] +
		-laavec[il3 + 1] * xx2[ix3]
	
	# only when k = 1         
	# inflow:
	# cladogenesis in one of the n-1 species: Q_M,n-1 -> Q_M,n; n+k-1 species present; rate once
	# extinction in one of the n+1 species: Q_M,n+1 -> Q_M,n; n+k+1 species present 
	# outflow:
	# all events with n+k species present
	dx3 = lacvec[il1] * nn[in4] * xx3[ix1] + muvec[il2] * nn[in2] * xx3[ix2] +
		-(lacvec[il3] + muvec[il3]) * nn[in3] * xx3[ix3] +
		-(laavec[il3] + gamvec[il3]) * xx3[ix3]
	
	return(list(c(dx1,dx2,dx3)))
}

library(odeintr)

integrate_sys(DAISIE_loglik_rhs())

# RHS1_odeintr --------------------------------------------------------------------


DAISIE_loglik_rhs_odeintr = function(t,x,pars)
{
  # Returns list with 3 strings containing right hand side of ODE for input in odeintr
  
  # Reads model parameters
  lx = (length(x) - 1)/2 # Number of variables containing the state of the system . x = current probs
  lac = pars[1]          # Lambda^c
  mu = pars[2]           # mu
  K = pars[3]            # K
  gam = pars[4]          # gamma
  laa = pars[5]          # Lambda^a
  kk = pars[6]           # k state (FILL HERE WHAT HIS MEANS)
  ddep = pars[7]         # Type of diversity dependence (0 = no DD; 1 = linear dep in speciation; 2 = exponential dep in 
                         # speciation rate; 11 = linear dep in speciation and immigration; 21 = exponential dep in speciation and immigration) 
  
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
  

  ############ COMPLETE THIS NEXT!


  listOfVectors <- list(laavec, lacvec, muvec, gamvec)
  listOfil <- list(il1, il2, il3, il4)
  vecOfVectorNames <- c("laavec", "lacvec", "muvec", "gamvec")
  `laavec[il1 + 1]` <- laavec[il1 + 1]
  `laavec[il3[1]]` <- laavec[il3[1]]
  `laavec[il3 + 1]` <- laavec[il3 + 1]
  `laavec[il3]` <- laavec[il3]
  `lacvec[il4 + 1]` <- lacvec[il4 + 1]
  `lacvec[il1]` <- lacvec[il1]
  `lacvec[il3]` <- lacvec[il3]
  `lacvec[il3[1]` <- lacvec[il3[1]]
  `lacvec[il1 + 1]` <- lacvec[il1 + 1]
  `lacvec[il3 + 1]` <- lacvec[il3 + 1]
  
  laavec[il1 + 1]
  
  # microbenchmark(
  # for(i in 1:length(vecOfVectorNames)){ 
  #   for(j in 1:length(listOfil)){
  #     assign(paste(vecOfVectorNames[[j]], "[il", i, "]", sep = "" ), listOfVectors[[i]][listOfil[[j]]])
  #   }
  # }
  # )
  
  temp_names <- names(temp_vector)
  paste("(", paste(temp_names, collapse = " + "), ")", sep = "")
  
  
  
  in1 = nil2lx+2*kk-1
  in2 = nil2lx+1
  in3 = nil2lx+kk
  
  ix1 = nil2lx-1
  ix2 = nil2lx+1
  ix3 = nil2lx
  ix4 = nil2lx-2

laavec[il1 + 1]  
laavec[il1 + 1] * xx2[ix1] + lacvec[il4 + 1] * xx2[ix4] + muvec[il2 + 1] * xx2[ix3] +
	lacvec[il1] * nn[in1] * xx1[ix1] + muvec[il2] * nn[in2] * xx1[ix2] +
	-(muvec[il3] + lacvec[il3]) * nn[in3] * xx1[ix3] +
	-gamvec[il3] * xx1[ix3]


# Do a paste function here that prints the sum of all the elements of each lac[il1]etc. vector.
# Output should be (lacvec[il1]1 + lacvec[il1]2)

# Attention: R does point to point vector multiplication. laa and xx2[ix1] are vectors
dx1_odeintr_char <- 'dxdt[0] = laa * x[1] + lambdaC * x[1] + mu * x[1] + lambdaC * n * x[0] + mu * n * x[0] - (mu + lambdaC) * n * x[0] - gam * x[0]; dxdt[1] = gam * x[0] + lambdaC * n * x[1] + mu * n * x[1] - (mu + lambdaC) * n * x[1] - laa * x[1]; dxdt[2] = -(laa + lambdaC + gam + mu) * x[2];'

# kk = 1
dx_odeintr_char <- 'dxdt[0] = laa * x[2] + 2 * lambdaC * x[2] + laa * x[1] + lambdaC * x[1] + mu * x[1] + lambdaC * n * x[0] + mu * n * x[0] - (mu + lambdaC) * n * x[0] - gam * x[0]; dxdt[1] = gam * x[0] + lambdaC * n * x[1] + mu * n * x[1] - (mu + lambdaC) * n * x[1] - laa * x[1]; dxdt[2] = lambdaC * n * x[2] + mu * n * x[2] - (lambdaC + mu) * n * x[2] - (laa + gam) * x[2];'

library(DAISIE)
library(odeintr)
pars = c(laa = laa, lambdaC = 0.1, mu = 0.1, gam = 0.1, n = 10)
pars = c(lacvec)

compile_sys("rhs1" , eq_system1, pars, TRUE)
compile_sys("rhs2" , eq_system2, pars, TRUE)

return(list(c(dx1,dx2,dx3)))
}


