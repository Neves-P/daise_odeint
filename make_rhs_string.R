#### Testing conditions ####
lac = .4
mu = 0.1
K = 30
gam = 0.2
laa = 0.3
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

##### Elements in  rhs ####

# laavec[il1+1]
# xx2[ix1]
# 
# lacvec[il4 + 1]
# xx2[ix4]
# 
# muvec[il2 + 1]
# xx2[ix3]
# 
# lacvec[il1]
# nn[in1]
# xx1[ix1]
# 
# muvec[il2]
# nn[in2]
# xx1[ix2]
# 
# -(muvec[il3] + lacvec[il3])
# nn[in3]
# xx1[ix3]
# 
# -gamvec[il3]
# xx1[ix3]

#### Function to make string ####

# for(i in 1:length(laavec)){
#   temp_str_1 <- paste("laa_il1_add_one", i, sep = "_")
#   
#   temp_str_2 <- paste("xx2_ix2", i, sep = "_")
#   laa_il1_add_one_xx2_ix1 <- append(laa_il1_add_one_xx2_ix1,
#                                     c(paste(temp_str_1, temp_str_2, sep = " * ")))
# }
# 
# for(i in 1:length(laavec[il1 + 1])){
#   temp_laavec <- paste("laavec_il1_plusone", i, sep = "_")
#   assign(temp_str1, laavec[il1 + 1][i])
#   temp_xx2 <- paste("xx2_ix1", i, sep = "_")
#   }
list_dx1 <- list()
dx1_pars_list <- list()
for(i in 1:length(laavec[il1])){
  temp_laavec_il1_plusone <- paste("laavec_il1_plusone", i, sep = "_")
  assign(temp_laavec_il1_plusone, laavec[il1 + 1][i])
  temp_xx2_ix1 <- paste("xx2_ix1", i, sep = "_")
  assign(temp_xx2_ix1, xx2[ix2][i])
  prod1 <- paste(temp_laavec_il1_plusone, temp_xx2_ix1, sep = " * ")
	print(prod1)
  temp_lacvec_il4_plusone <- paste("lacvec_il4_plusone", i, sep = "_")
  assign(temp_lacvec_il4_plusone, lacvec[il4 + 1][i])
  temp_xx2_ix4 <- paste("xx2_ix4", i, sep = "_")
  assign(temp_xx2_ix4, xx2[ix4][i])
  prod2 <- paste(temp_lacvec_il4_plusone, temp_xx2_ix4, sep = " * ")
  
  temp_muvec_il2_plusone <- paste("muvec_il2_plusone", i, sep = "_")
  assign(temp_muvec_il2_plusone, muvec[il2 + 1][i])
  temp_xx2_ix3 <- paste("xx2_ix3", i, sep = "_")
  assign(temp_xx2_ix3, xx2[ix3][i])
  prod3 <- paste(temp_muvec_il2_plusone, temp_xx2_ix3, sep = " * ")
  
  temp_lacvec_il1 <- paste("lacvec_il1", i, sep = "_")
  assign(temp_lacvec_il1, lacvec[il1][i])
  temp_nn_in1 <- paste("nn_in1", i, sep = "_")
  assign(temp_nn_in1, nn[in1][i])
  temp_xx1_ix1 <- paste("xx1_ix1", i, sep = "_")
  assign(temp_xx1_ix1, xx1[ix1][i])
  prod4 <- paste(temp_lacvec_il1, temp_nn_in1, temp_xx1_ix1, sep = " * ")
  
  temp_muvec_il3 <- paste("muvec_il3", i, sep= "_")
  assign(temp_muvec_il3, muvec[il3][i])
  temp_lacvec_il3 <- paste("lacvec_il3", i, sep = "_")
  assign(temp_lacvec_il3, lacvec[il3][i])
  neg_term1 <- paste("-(", paste(temp_muvec_il3, temp_lacvec_il3, 
                                 sep = " + "), ")", sep = "")

  temp_nn_in2 <- paste("nn_in2", i, sep = "_")
  assign(temp_nn_in2, nn[in2][i])
  temp_xx1_ix3 <- paste("xx1_ix3", i, sep = "_")
  assign(temp_xx1_ix3, xx1[ix3][i])
  prod5 <- paste(neg_term1, temp_nn_in2, temp_xx1_ix3, sep = " * ")
  
  temp_neggamvec_il3 <- paste("-gamvec_il3", i, sep = "_")
  assign(temp_neggamvec_il3, -gamvec[il3][i]) # GAMVEC CHANGED SIGN
  prod6 <- paste(temp_neggamvec_il3, temp_xx1_ix3, sep = " * ")
  
  complete_rhs <- paste(prod1, prod2, prod3, prod4, prod5, prod6, sep = " + ")
  list_dx1[[i]] <- complete_rhs

  
  #Pars

	pars    <- ls(pattern = "temp")
	pars <- pars[grepl("temp", pars) & !grepl("xx", pars)]
  dx1_pars_list[[i]] <- mget(unlist(unname(mget(pars)))) # Recovers parameters into list
}



 