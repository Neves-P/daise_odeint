#### Testing conditions #### CHECK IF LIST_INDICES WORKS
lac = pars[1]
mu = pars1[2]
K = pars1[3]
gam = pars1[4]
laa = pars1[5]
kk = 10
ddep = 0
x <- c(1,2,3,4,5,6,7,8,9,10)
lx = (length(x) - 1)/2

nn = -2:(lx+2*kk+1)
lnn = length(nn)         # 
nn = pmax(rep(0,lnn),nn) # Number of species in clade

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
	
	list_dx <- list()
	pars_list <- list()
	par_name_list <- list()
	init_state_list <- list()
	init_state_list_names <-list()
	x_counter <- 0
	dx_list_counter <- 1
	
	for(i in 1:length(list_pars$laavec[list_indices$il1])){
		
		#### dx1 ####
		
		# First product
		temp_laavec_il1_plusone <- paste("laavec_il1_plusone", i, sep = "_")
		assign(temp_laavec_il1_plusone, list_pars$laavec[list_indices$il1+1][i])

		temp_xx2_ix1 <- paste0("x[", x_counter, "]")
		assign(temp_xx2_ix1, list_pars$xx2[list_indices$ix1][i])
		
		prod1 <- paste(temp_laavec_il1_plusone, temp_xx2_ix1, sep = " * ")
		
		
		# Second product
		temp_lacvec_il4_plusone <- paste("lacvec_il4_plusone", i, sep = "_")
		assign(temp_lacvec_il4_plusone, lacvec[il4 + 1][i])
		
		temp_xx2_ix4 <- paste0("x[", x_counter + 1, "]")
		assign(temp_xx2_ix4, xx2[ix4][i])
		
		prod2 <- paste(temp_lacvec_il4_plusone, temp_xx2_ix4, sep = " * ")
		
		
		# Third product
		temp_muvec_il2_plusone <- paste("muvec_il2_plusone", i, sep = "_")
		assign(temp_muvec_il2_plusone, list_pars$muvec[list_indices$il2 + 1][i])
		
		temp_xx2_ix3 <- paste0("x[", x_counter + 2, "]")
		assign(temp_xx2_ix3, list_pars$xx2[list_indices$ix3][i])
		
		prod3 <- paste(temp_muvec_il2_plusone, temp_xx2_ix3, sep = " * ")

		
		# Fourth product
		temp_lacvec_il1 <- paste("lacvec_il1", i, sep = "_")
		assign(temp_lacvec_il1, list_pars$lacvec[list_indices$il1][i])
		
		temp_nn_in1 <- paste("nn_in1", i, sep = "_")
		assign(temp_nn_in1, list_pars$nn[list_indices$in1][i])
		
		temp_xx1_ix1 <- paste0("x[", x_counter + 3, "]")
		assign(temp_xx1_ix1, list_pars$xx1[list_indices$ix1][i])
		
		prod4 <- paste(temp_lacvec_il1, temp_nn_in1, temp_xx1_ix1, sep = " * ")
		
		
		
		# Actual fifth product
		temp_muvec_il2 <- paste("muvec_il2", i, sep = "_")
		assign(temp_muvec_il2, list_pars$muvec[list_indices$il2][i])
		
		temp_nn_in2 <- paste("nn_in2", i, sep = "_")
		assign(temp_nn_in2, list_pars$nn[list_indices$in2][i])
		
		temp_xx1_ix2 <- paste0("x[", x_counter + 4, "]")
		assign(temp_xx1_ix2, list_pars$xx1[list_indices$ix2][i])
		
		prod5 <- paste(temp_muvec_il2, temp_nn_in2, temp_xx1_ix2, sep = " * ")
		
		
		
		# Negative term
		temp_muvec_il3 <- paste("muvec_il3", i, sep= "_")
		assign(temp_muvec_il3, list_pars$muvec[list_indices$il3][i])
		
		temp_lacvec_il3 <- paste("lacvec_il3", i, sep = "_")
		assign(temp_lacvec_il3, list_pars$lacvec[list_indices$il3][i])
		
		neg_term1 <- paste("-(", paste(temp_muvec_il3, temp_lacvec_il3, 
																	 sep = " + "), ")", sep = "")
		
		
		# Sixth product
		temp_nn_in2 <- paste("nn_in2", i, sep = "_")
		assign(temp_nn_in2, list_pars$nn[list_indices$in2][i])
		
		temp_xx1_ix3 <- paste0("x[", x_counter + 5, "]")
		assign(temp_xx1_ix3, list_pars$xx1[list_indices$ix3][i])
		
		prod6 <- paste(neg_term1, temp_nn_in2, temp_xx1_ix3, sep = " * ")
		
		
		# Seventh product
		# This gamvec is assign the regular (non negative) value of gam
		# The string is built with a negative sign for correct input in odeintr
		temp_neggamvec_il3 <- paste("- gamvec_il3", i, sep = "_")
		assign(temp_neggamvec_il3, list_pars$gamvec[list_indices$il3][i]) 
		
		prod7 <- paste(temp_neggamvec_il3, temp_xx1_ix3, sep = " * ")
		
		# dx1 rhs of equation
		complete_rhs <- paste(prod1, prod2, prod3, prod4, prod5, prod6, prod7, sep = " + ")
		list_dx[[dx_list_counter]] <- complete_rhs
		
		
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
		
		temp_xx2_ix2 <- paste0("x[", x_counter + 2, "]")
		assign(temp_xx2_ix2, list_pars$xx2[list_indices$ix2][i])
		
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
		
		
		prod5 <- paste(temp_laavec_il3_plusone, temp_xx2_ix3, sep = " * ")
		
		# dx2 rhs of equation
		complete_rhs <- paste(prod1, prod2, prod3, prod4, prod5, sep = " + ")
		list_dx[[dx_list_counter + 1]] <- complete_rhs
		
		
		#### dx3 ####
		
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
	
		prod1 <- paste(neg_term1, paste0("x[", x_counter + 6, "]"), sep = " * ") 
		
		list_dx[[dx_list_counter + 2]] <- prod1
		
	
		#### Model parameters per rhs ####
		
		# Updates index of X and dx_list for next equation loop
		
		x_counter <- x_counter + 7
		dx_list_counter <- dx_list_counter + 3
		
		# Store initial state and parameters
		local_env_pars <- ls(pattern = "temp")
		
		pars <- local_env_pars[grepl("temp", local_env_pars) &
													 	!grepl("xx", local_env_pars)]
		par_name_list[[i]] <- unlist(unname(mget(pars)))
		pars_list[[i]] <- mget(unlist(unname(mget(pars)))) # Recovers parameters into list
		
		init_state_list_names[[i]] <- mget(local_env_pars[grepl("temp", local_env_pars) 
																 & grepl("xx", local_env_pars)])
		init_state_list[[i]] <- mget(unlist(unname(init_state_list_names[[i]])))
	}
	# Return
	
	return(list(rhs = list_dx, pars = unlist(pars_list), 
							init_state = unlist(init_state_list)))
}

make_rhs_1(list_pars = list_pars, list_indices = list_indices)

make_sys <- function(rhs)
	{
	ode_system <- list()
	for(i in 1:length(rhs$rhs)){
		ode_system[[i]] <- (paste0("dxdt[", i - 1, "] = ", rhs$rhs[[i]], "; "))
	}

	return(paste(ode_system, sep = "", collapse = ""))
}

sys <- make_rhs_1(list_pars, list_indices)
eqs <- make_sys(sys)

write(tail(eqs), "dxdt.txt")
y <- compile_sys(name = "y", sys$rhs, pars[unique(names(pars))], sys_dim = length(x)) # Elements of

write(y, "the_code_new.cpp")
pars_list_names <- make_rhs_1(list_pars, list_indices)
pars <- pars_list_names$pars
pars[unique(names(pars))]


compile_sys(name = "y", sys$rhs, pars[unique(names(pars))])
list_pars


compiled <- y(x, 10, 0.01)
initial <- c(rep(0, 9))
length(initial)
the_code <- y
write(the_code, "the_code")


