# Help ####
# Basics
# http://cran.r-project.org/doc/contrib/Leisch-CreatingPackages.pdf

# Building package skeleton ####

# List all the code files
code_files <-list.files("./R", full.names=TRUE)

# Build package skeleton
package.skeleton(name="dendro-toolkit", code_files=code_files, force=TRUE)