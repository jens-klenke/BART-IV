# load packages 
############## Packages ################
source(here::here('01_code/packages.R'))
# some packages have to be implemented by hand (from github)

# source all files in the functions folder
invisible(sapply(list.files(here::here('01_code/functions'), full.names = TRUE), 
                 source))

## load data
### Data to big for Github -> download from -> Sciebo/BART-IV-Data 
load(url('https://uni-duisburg-essen.sciebo.de/s/1kkmrWzX38EL7ID/download'))





