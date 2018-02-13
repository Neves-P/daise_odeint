# Writen by Pedro Neves on 15/12/17, under the GPL-2 license.
# Code adapted from package DAISIE (Etienne, Valente, Phillimore & Haegeman), 
# requires package odeintr (Keitt)

# TODO(Neves-P): Complete compile_sys or integrate_sys input function

require(odeintr)

DAISIE_loglik = function(pars1,pars2,brts,stac,missnumspec,methode = "lsodes", solver = "deSolve")
{
  # Allows the use of two different ODE integrators
  
  # brts = branching times (positive, from present to past)
  # - max(brts) = age of the island
  # - next largest brts = stem age / time of divergence from the mainland
  # The interpretation of this depends on stac (see below)
  # For stac = 0, there is no other value.
  # For stac = 1 and stac = 5, this is the time since divergence from the immigrant's sister on the mainland.
  # The immigrant must have immigrated at some point since then.
  # For stac = 2 and stac = 3, this is the time since divergence from the mainland.
  # The immigrant that established the clade on the island must have immigrated precisely at this point.
  # For stac = 3, it must have reimmigrated, but only after the first immigrant had undergone speciation.
  # - min(brts) = most recent branching time (only for stac = 2, or stac = 3)
  # pars1 = model parameters
  # - pars1[1] = lac = (initial) cladogenesis rate
  # - pars1[2] = mu = extinction rate
  # - pars1[3] = K = maximum number of species possible in the clade
  # - pars1[4] = gam = (initial) immigration rate
  # - pars1[5] = laa = (initial) anagenesis rate
  # pars2 = model settings
  # - pars2[1] = lx = length of ODE variable x
  # - pars2[2] = ddep = diversity-dependent model,mode of diversity-dependence
  #  . ddep == 0 : no diversity-dependence
  #  . ddep == 1 : linear dependence in speciation rate (anagenesis and cladogenesis)
  #  . ddep == 11 : linear dependence in speciation rate and immigration rate
  #  . ddep == 3 : linear dependence in extinction rate
  # - pars2[3] = cond = conditioning
  #  . cond == 0 : no conditioning
  #  . cond == 1 : conditioning on presence on the island (not used in this single loglikelihood)
  # - pars2[4] = parameters and likelihood should be printed (1) or not (0)
  # stac = status of the clade formed by the immigrant
  #  . stac == 0 : immigrant is not present and has not formed an extant clade
  #  . stac == 1 : immigrant is present but has not formed an extant clade
  #  . stac == 2 : immigrant is not present but has formed an extant clade
  #  . stac == 3 : immigrant is present and has formed an extant clade
  #  . stac == 4 : immigrant is present but has not formed an extant clade, and it is known when it immigrated.
  #  . stac == 5 : immigrant is not present and has not formed an extant clade, but only an endemic species
  # missnumspec = number of missing species
  # odeinter == FALSE : Use default ODE solver deSolve. If TRUE use odeintr
  
  ddep = pars2[2]
  cond = pars2[3]
  if(cond > 0)
  {
    cat("Conditioning has not been implemented and may not make sense. Cond is set to 0.\n")
  }
  
  lac = pars1[1]
  mu = pars1[2]
  K = pars1[3]
  if(ddep == 0)
  {
    K = Inf
  }
  gam = pars1[4]
  laa = pars1[5]
  
  abstol = 1e-16
  reltol = 1e-10
  brts = -sort(abs(as.numeric(brts)),decreasing = TRUE)
  if(length(brts) == 1 & sum(brts == 0) == 1)
  {
    stop('The branching times contain only a 0. This means the island emerged at the present which is not allowed.');
    loglik = -Inf
    return(loglik)
  }
  if(sum(brts == 0) == 0)
  {
    brts[length(brts) + 1] = 0
  }
  
  if(solver != "deSolve" || solver != "odeintr"){
    stop('Specified ODE solver not recognised. Please select "deSolve" or "odeintr"')
  }
  
  if(solver == "deSolve"){
    
    # for stac = 0 and stac = 1, brts will contain origin of island and 0; length = 2; no. species should be 0
    # for stac = 1, brts will contain origin of island and 0; length = 2; no. species should be 1
    # for stac = 4, brts will contain origin of island, colonization event and 0; length = 3; no. species should be 1
    # for stac = 2, brts with contain origin of island, colonization event, branching times, 0; no. species should be no. branching times + 1
    # for stac = 3, brts with contain origin of island, colonization event, branching times, 0; no. species should be no. branching times + 2
    #S = length(brts) - 2 * (stac == 0) - (stac == 1) - 2 * (stac == 4) - 2 * (stac == 2) - (stac == 3)
    S = length(brts) - (stac %% 2 == 1) - 2 * (stac %% 2 == 0)
    S2 = S - (stac == 1) - (stac == 3) - (stac == 4)
    loglik = -lgamma(S2 + missnumspec + 1) + lgamma(S2 + 1) + lgamma(missnumspec + 1)
    if(min(pars1) < 0)
    {
      cat('One or more parameters are negative.\n')
      loglik = -Inf
      return(loglik)
    }
    if((ddep == 1 | ddep == 11) & ceiling(K) < (S + missnumspec))
    {
      #cat('The proposed value of K is inompatible with the number of species in the clade. Likelihood for this parameter set will be set to -Inf.\n')
      loglik = -Inf
      return(loglik)
    }
    if(lac == Inf & mu != Inf & missnumspec == 0)
    {
      loglik = DAISIE_loglik_high_lambda(pars1,-brts,stac)
    } else {
      if(ddep == 1 | ddep == 11)
      {
        lx = min(1 + max(missnumspec,ceiling(K)),round(pars2[1]) + missnumspec)
      } else {
        lx = roundn(pars2[1]) + missnumspec
      }
      if(loglik > -Inf)
      { 
        # in all cases we integrate from the origin of the island to the first branching point (stac > 1) or to the present (stac <= 1)
        probs = rep(0,2 * lx + 1)
        probs[1] = 1
        k1 = 0
        y = ode(probs,brts[1:2],DAISIE_loglik_rhs,c(pars1,k1,ddep),rtol = reltol,atol = abstol,method = methode)
        probs = y[2,2:(2 * lx + 2)]
        cp = checkprobs(lv = 2 * lx,loglik,probs); loglik = cp[[1]]; probs = cp[[2]]      
        if(stac == 0)
          # for stac = 0, the integration is from the origin of the island until the present
          # and we evaluate the probability of no clade being present and no immigrant species,
          # but there can be missing species
        {     
          loglik = loglik + log(probs[1 + missnumspec])
        } else {
          if(stac == 1 || stac == 5)
            # for stac = 1, the integration is from the maximum colonization time (usually the
            # island age + tiny time unit) until the present, where we set all probabilities where
            # the immigrant is already present to 0
            # and we evaluate the probability of the immigrant species being present,
            # but there can be missing species 
            # for stac = 5, we do exactly the same, but we evaluate the probability of an endemic species being present alone.          
          {         
            probs[(lx + 1):(2 * lx)] = 0
            y = ode(probs,brts[2:3],DAISIE_loglik_rhs,c(pars1,k1,ddep),rtol = reltol,atol = abstol,method = methode)
            probs = y[2,2:(2 * lx + 2)]
            cp = checkprobs(lv = 2 * lx,loglik,probs); loglik = cp[[1]]; probs = cp[[2]]               
            loglik = loglik + log(probs[(stac == 1) * lx + (stac == 5) + 1 + missnumspec])
          } else {
            # for stac > 1, but not 5, integration is then from the colonization event until the first branching time (stac = 2 and 3) or the present (stac = 4). We add a set of equations for Q_M,n, the probability that the process is compatible with the data, and speciation has not happened; during this time immigration is not allowed because it would alter the colonization time. After speciation, colonization is allowed again (re-immigration)
            # all probabilities of states with the immigrant present are set to zero and all probabilities of states with endemics present are transported to the state with the colonist present waiting for speciation to happen. We also multiply by the (possibly diversity-dependent) immigration rate
            gamvec = divdepvec(gam,K,lx,k1,ddep * (ddep == 11 | ddep == 21))
            probs[(2 * lx + 1):(3 * lx)] = gamvec[1:lx] * probs[1:lx]
            probs[1:(2 * lx)] = 0        
            k1 = 1
            y = ode(probs,c(brts[2:3]),DAISIE_loglik_rhs2,c(pars1,k1,ddep),rtol = reltol,atol = abstol,method = methode)
            probs = y[2,2:(3 * lx + 1)]
            cp = checkprobs2(lx,loglik,probs); loglik = cp[[1]]; probs = cp[[2]]
            if(stac == 4)
              # if stac = 4, we're done and we take an element from Q_M,n
            {
              loglik = loglik + log(probs[2 * lx + 1 + missnumspec])
            } else {         
              # for stac = 2 and 3, at the first branching point all probabilities of states Q_M,n are transferred to probabilities where only endemics are present. Then go through the branching points.
              S1 = length(brts) - 1
              if(S1 >= 3)
              {
                lacvec = divdepvec(lac,K,lx,k1,ddep)
                probs[1:lx] = lacvec[1:lx] * (probs[1:lx] + probs[(2 * lx + 1):(3 * lx)])
                probs[(lx + 1):(2 * lx)] = lacvec[2:(lx + 1)] * probs[(lx + 1):(2 * lx)]
                probs = probs[-c((2 * lx + 2):(3 * lx))]
                probs[2 * lx + 1] = 0
                for(k in 3:S1)
                {
                  k1 = k - 1
                  y = ode(probs,brts[k:(k+1)],DAISIE_loglik_rhs,c(pars1,k1,ddep),rtol = reltol,atol = abstol,method = methode)
                  probs = y[2,2:(2 * lx + 2)]
                  if(k < S1)
                  {
                    # speciation event      
                    lacvec = divdepvec(lac,K,lx,k1,ddep)
                    probs[1:(2 * lx)] = c(lacvec[1:lx],lacvec[2:(lx + 1)]) * probs[1:(2 * lx)]
                  }
                }            
              }
              # we evaluate the probability of the phylogeny with any missing species at the present without (stac = 2) or with (stac = 3) the immigrant species; there can be no missing species for stac = 4
              loglik = loglik + log(probs[(stac == 3) * lx + 1 + missnumspec])
            }   
          }     
        }           
      }
    }
  }
  
  
  if(solver == "odeintr"){
    # for stac = 0 and stac = 1, brts will contain origin of island and 0; length = 2; no. species should be 0
    # for stac = 1, brts will contain origin of island and 0; length = 2; no. species should be 1
    # for stac = 4, brts will contain origin of island, colonization event and 0; length = 3; no. species should be 1
    # for stac = 2, brts with contain origin of island, colonization event, branching times, 0; no. species should be no. branching times + 1
    # for stac = 3, brts with contain origin of island, colonization event, branching times, 0; no. species should be no. branching times + 2
    #S = length(brts) - 2 * (stac == 0) - (stac == 1) - 2 * (stac == 4) - 2 * (stac == 2) - (stac == 3)
    S = length(brts) - (stac %% 2 == 1) - 2 * (stac %% 2 == 0)
    S2 = S - (stac == 1) - (stac == 3) - (stac == 4)
    loglik = -lgamma(S2 + missnumspec + 1) + lgamma(S2 + 1) + lgamma(missnumspec + 1)
    if(min(pars1) < 0)
    {
      cat('One or more parameters are negative.\n')
      loglik = -Inf
      return(loglik)
    }
    if((ddep == 1 | ddep == 11) & ceiling(K) < (S + missnumspec))
    {
      #cat('The proposed value of K is incompatible with the number of species in the clade. Likelihood for this parameter set will be set to -Inf.\n')
      loglik = -Inf
      return(loglik)
    }
    if(lac == Inf & mu != Inf & missnumspec == 0)
    {
      loglik = DAISIE_loglik_high_lambda(pars1,-brts,stac)
    } else {
      if(ddep == 1 | ddep == 11)
      {
        lx = min(1 + max(missnumspec,ceiling(K)),round(pars2[1]) + missnumspec)
      } else {
        lx = roundn(pars2[1]) + missnumspec
      }
      if(loglik > -Inf)
      { 
        # in all cases we integrate from the origin of the island to the first branching point (stac > 1) or to the present (stac <= 1)
        probs = rep(0, 2 * lx + 1)
        probs[1] = 1
        k1 = 0
        y = compile_sys("y", DAISIE_loglik_rhs_odeintr(probs, c(pars1, k1, ddep)))
        y = ode(probs,brts[1:2],DAISIE_loglik_rhs,c(pars1,k1,ddep),rtol = reltol,atol = abstol,method = methode)
        probs = y[2,2:(2 * lx + 2)]
        cp = checkprobs(lv = 2 * lx,loglik,probs); loglik = cp[[1]]; probs = cp[[2]]      
        if(stac == 0)
          # for stac = 0, the integration is from the origin of the island until the present
          # and we evaluate the probability of no clade being present and no immigrant species,
          # but there can be missing species
        {     
          loglik = loglik + log(probs[1 + missnumspec])
        } else {
          if(stac == 1 || stac == 5)
            # for stac = 1, the integration is from the maximum colonization time (usually the
            # island age + tiny time unit) until the present, where we set all probabilities where
            # the immigrant is already present to 0
            # and we evaluate the probability of the immigrant species being present,
            # but there can be missing species 
            # for stac = 5, we do exactly the same, but we evaluate the probability of an endemic species being present alone.          
          {         
            probs[(lx + 1):(2 * lx)] = 0
            y = ode(probs,brts[2:3],DAISIE_loglik_rhs,c(pars1,k1,ddep),rtol = reltol,atol = abstol,method = methode)
            probs = y[2,2:(2 * lx + 2)]
            cp = checkprobs(lv = 2 * lx,loglik,probs); loglik = cp[[1]]; probs = cp[[2]]               
            loglik = loglik + log(probs[(stac == 1) * lx + (stac == 5) + 1 + missnumspec])
          } else {
            # for stac > 1, but not 5, integration is then from the colonization event until the first branching time (stac = 2 and 3) or the present (stac = 4). We add a set of equations for Q_M,n, the probability that the process is compatible with the data, and speciation has not happened; during this time immigration is not allowed because it would alter the colonization time. After speciation, colonization is allowed again (re-immigration)
            # all probabilities of states with the immigrant present are set to zero and all probabilities of states with endemics present are transported to the state with the colonist present waiting for speciation to happen. We also multiply by the (possibly diversity-dependent) immigration rate
            gamvec = divdepvec(gam,K,lx,k1,ddep * (ddep == 11 | ddep == 21))
            probs[(2 * lx + 1):(3 * lx)] = gamvec[1:lx] * probs[1:lx]
            probs[1:(2 * lx)] = 0        
            k1 = 1
            y = ode(probs,c(brts[2:3]),DAISIE_loglik_rhs2,c(pars1,k1,ddep),rtol = reltol,atol = abstol,method = methode)
            probs = y[2,2:(3 * lx + 1)]
            cp = checkprobs2(lx,loglik,probs); loglik = cp[[1]]; probs = cp[[2]]
            if(stac == 4)
              # if stac = 4, we're done and we take an element from Q_M,n
            {
              loglik = loglik + log(probs[2 * lx + 1 + missnumspec])
            } else {         
              # for stac = 2 and 3, at the first branching point all probabilities of states Q_M,n are transferred to probabilities where only endemics are present. Then go through the branching points.
              S1 = length(brts) - 1
              if(S1 >= 3)
              {
                lacvec = divdepvec(lac,K,lx,k1,ddep)
                probs[1:lx] = lacvec[1:lx] * (probs[1:lx] + probs[(2 * lx + 1):(3 * lx)])
                probs[(lx + 1):(2 * lx)] = lacvec[2:(lx + 1)] * probs[(lx + 1):(2 * lx)]
                probs = probs[-c((2 * lx + 2):(3 * lx))]
                probs[2 * lx + 1] = 0
                for(k in 3:S1)
                {
                  k1 = k - 1
                  y = ode(probs,brts[k:(k+1)],DAISIE_loglik_rhs,c(pars1,k1,ddep),rtol = reltol,atol = abstol,method = methode)
                  probs = y[2,2:(2 * lx + 2)]
                  if(k < S1)
                  {
                    # speciation event      
                    lacvec = divdepvec(lac,K,lx,k1,ddep)
                    probs[1:(2 * lx)] = c(lacvec[1:lx],lacvec[2:(lx + 1)]) * probs[1:(2 * lx)]
                  }
                }            
              }
              # we evaluate the probability of the phylogeny with any missing species at the present without (stac = 2) or with (stac = 3) the immigrant species; there can be no missing species for stac = 4
              loglik = loglik + log(probs[(stac == 3) * lx + 1 + missnumspec])
            }   
          }     
        }           
      }
    }
  } 
}
