dll_folders <- list.dirs(paste0(tempdir(), "\\sourceCpp-x86_64-w64-mingw32-0.12.15"))

dll_folders <- dll_folders[2:length(dll_folders)]

paste0(dll_folders[[1]], "\\", dir(dll_folders[[1]]))

full_paths <- c()
for(i in 1:length(dll_folders)){
  full_paths[i] <- paste0(dll_folders[[i]], "\\", dir(dll_folders[[i]]))
}

loaded_dll <- getLoadedDLLs()

loaded_dll_cpp <- loaded_dll[grep("sourceCpp", names_dll_list)]

paths <- c()
for (i in 1:length(loaded_dll_cpp)) {
  paths[i] <- loaded_dll_cpp[i][[1]][[2]]  
  
}
for (i in 5:length(paths)){
  dyn.unload(paths[[i]])
}
getLoadedDLLs()
