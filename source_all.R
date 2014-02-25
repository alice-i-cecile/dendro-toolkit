# Source all package files
R_files <- list.files("./R", full.names=T)

lapply(R_files, source)
