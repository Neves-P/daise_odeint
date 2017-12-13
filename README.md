# daise_odeint

[![Join the chat at https://gitter.im/daise_odeint/Lobby](https://badges.gitter.im/daise_odeint/Lobby.svg)](https://gitter.im/daise_odeint/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)



daisie_odeint is a simple rewrite of the daise_loglik_rhs and daisie_loglik_rhs2 internal functions from the R package DAISIE: Dynamical Assembly of Islands by Speciation, Immigration and Extinction, by Etienne, Valente, Phillimore, and Haegeman. This pair of functions is called by the main function chain (daisie_ML, daisie_loglik_all).

daisie_odeint allows the use of another integration engine, the C++ library odeint from the Boost collection, apart from the already implemented [deSolve](https://cran.r-project.org/package=deSolve). 
To achieve this, daisie_odeint relies on the R package [odeintr](https://cran.r-project.org/package=odeintr), a wrapper for odeint in R.
Possible custom Rcpp integration can be implemented in the future

The original version of DAISIE can be found on CRAN, [here](https://cran.r-project.org/package=DAISIE).
