# daise_odeint


daisie_odeint is a simple rewrite of the daise_loglik_rhs and daisie_loglik_rhs2 internal functions from the R package DAISIE: Dynamical Assembly of Islands by Speciation, Immigration and Extinction, by Etienne, Valente, Phillimore, and Haegeman. This pair of functions is called by the main function chain (daisie_ML -> daisie_loglik_all).

daisie_odeint allows the use of another integration engine, the C++ library odeint from the Boost collection, apart from the already implemented [deSolve](https://cran.r-project.org/package=deSolve). 
To achieve this, daisie_odeint relies on the R package [odeintr](https://cran.r-project.org/package=odeintr), a wrapper for odeint in R.
Possible custom Rcpp integration can be implemented in the future.


The original version of DAISIE can be found on CRAN, [here](https://cran.r-project.org/package=DAISIE).

daisie_odeint is not yet correctly compiling and working with odeintr.
See a roadmap of the project [here](https://github.com/Neves-P/daisie_odeint/projects/1).


# References
1. Valente, Luis M., Albert B. Phillimore, and Rampal S. Etienne. "Equilibrium and non‐equilibrium dynamics simultaneously operate in the Galápagos islands." Ecology letters 18.8 (2015): 844-852. [doi:10.1111/ele.12461](https://doi.org/10.1111/ele.12461)

2. Soetaert, Karline, Thomas Petzoldt, and R Woodrow Setzer. 2010. “Solving Differential Equations in R : Package deSolve.” Journal of Statistical Software 33 (9): 1–25. doi:10.18637/jss.v033.i09.

3. Timothy H. Keitt (2017). odeintr: C++ ODE Solvers Compiled on-Demand. R package version 1.7.1. https://CRAN.R-project.org/package=odeintr

4. Ahnert, Karsten, and Mario Mulansky. 2011. “Odeint – Solving Ordinary Differential Equations in C++.” In AIP Conference Proceedings, 1389:1586–89. American Institute of Physics. doi:10.1063/1.3637934.

