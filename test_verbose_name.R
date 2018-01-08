xx1
prod <- lacvec[il4 + 1] * xx2[ix4]


for(i in 1:length(prod)){
  nam <- paste("prod", i, sep = "")
  assign(nam, prod[i])
}


for(i in 1:length(lacvec[il4 + 1])){
  nam <- paste("lacvec[il4 + 1]", i, sep = "")
  assign(nam, lacvec[il4 + 1][i])
}

for(i in 1:length(xx2[ix4])){
  nam <- paste("xx2[ix4]", i, sep = "")
  assign(nam, xx2[ix4][i])
}

lac_il4_add_one_xx2_ix4 <- c()
for(i in 1:length(lacvec[il4 + 1])){
  temp_str_1 <- paste("lacvec[il4 + 1]", i, sep = "_")
  temp_str_2 <- paste("xx2[ix4]", i, sep = "_")
  lac_il4_add_one_xx2_ix4 <- append(lac_il4_add_one_xx2_ix4, 
                                    c(paste(temp_str_1, temp_str_2, sep = " * ")))
}


#laa[il1 + 1] * xx2[ix1] string segment *Vectors are not the same size!*
laa_il1_add_one_xx2_ix1 <- c()
for(i in 1:length(laavec)){
	temp_str_1 <- paste("laa_il1_add_one", i, sep = "_")
	
	temp_str_2 <- paste("xx2_ix2", i, sep = "_")
	laa_il1_add_one_xx2_ix1 <- append(laa_il1_add_one_xx2_ix1, c(paste(temp_str_1, temp_str_2, sep = " * ")))
}




