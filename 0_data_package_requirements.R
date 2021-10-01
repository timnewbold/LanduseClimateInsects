##%######################################################%##
#                                                          #
####               R package requirements               ####
#                                                          #
##%######################################################%##

# This script details all of the datasets and R packages required to 
# run the analyses associated with the paper by Outhwaite et al on the
# interactions between land use/use intensity and climate change on insects.

############ 1. R packages required from github ############

# A number of R packages have been developed by Dr Tim Newbold for use when
# analysing the PREDICTS database. Code to download and install these 
# packages from Github is detailed here.

# the devtools library is required to download packages from Github
library(devtools)

# The library predictsFunctions includes various functions for organising, analysing
# and plotting outputs from analyses of the PREDICTS database. 
install_github(repo = "timnewbold/predicts-demo",subdir = "predictsFunctions")

# The StatisticalModels package includes various functions for analysing the
# PREDICTS database
install_github(repo = "timnewbold/StatisticalModels")

