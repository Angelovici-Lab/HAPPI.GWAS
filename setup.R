## Clean the workspace
rm(list = ls())



# Set repository so that required packages can be downloaded
r = getOption("repos")
r["CRAN"] = "https://cloud.r-project.org/"
options(repos = r)



R_library = paste(R.version$platform, "-library", sep = "")
R_version = gsub("\\.0$", "", paste(R.version$major, R.version$minor, sep = "."))

# Install path
p <- file.path("~/R", R_library, R_version)
if(!dir.exists(p)){
  dir.create(path = p, recursive = TRUE)
}



# Gather required packages
packages <- c(
    "ape",
    "BiocManager",
    "compiler", "car",
    "data.table", "DataCombine",
    "EMMREML",
    "genetics", "gplots", "gridExtra",
    "lme4", "LDheatmap",
    "scatterplot3d",
    "dplyr", "tidyr", "tibble", "ggplot2",
    "yaml",
    "foreach", "doParallel",
    "bigmemory", "biganalytics"
)

# Check packages and install them if needed
invisible(lapply(packages, FUN = function(x){
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE, repos = "https://cloud.r-project.org/", lib = p)
    library(x, lib.loc = p, character.only = TRUE)
  }
}))



# The packages here are from BiocManager
bioc_packages <- c("Biobase", "BiocGenerics", "zlibbioc", "snpStats", "multtest")

# Check packages and install them if needed
invisible(lapply(bioc_packages, FUN = function(x){
  if (!require(x, character.only = TRUE)) {
    BiocManager::install(x, lib.loc = p, lib = p, update = TRUE, ask = FALSE)
    library(x, lib.loc = p, character.only = TRUE)
  }
}))



# Print this after all packages are successfully installed
loaded_packages <- sessionInfo()
print(names(loaded_packages$otherPkgs))



# Print some empty space
cat(rep("\n", 3))



# Create these folders
foldernames <- c("raw_data", "BLUP_BLUE", "reference_files", "output", "yaml")
for (i in 1:length(foldernames)) {
  temp <- file.path("..", foldernames[i])
  if(!dir.exists(temp)){
    dir.create(temp)
    if(dir.exists(temp)){
      print(paste(foldernames[i], " folder has been created!!!", sep = ""))
    } else{
      print(paste(foldernames[i], " folder cannot be created!!!", sep = ""))
    }
  } else{
    print(paste(foldernames[i], " folder exists!!!", sep = ""))
  }
}



# Create Arabidopsis360.yaml file
filename <- "Arabidopsis360.yaml"
dat <- readLines(con = filename)
writeLines(dat, con = file.path("..", filename))
if (file.exists(file.path("..", filename))) {
  print(paste(filename, " has been created!!!", sep = ""))
} else{
  print(paste(filename, " cannot be created!!!", sep = ""))
}



# Create Arabidopsis1001.yaml file 
filename <- "Arabidopsis1001.yaml"
dat <- readLines(con = filename)
writeLines(dat, con = file.path("..", filename))
if (file.exists(file.path("..", filename))) {
  print(paste(filename, " has been created!!!", sep = ""))
} else{
  print(paste(filename, " cannot be created!!!", sep = ""))
}



# Create Maize.yaml file
filename <- "Maize.yaml"
dat <- readLines(con = filename)
writeLines(dat, con = file.path("..", filename))
if (file.exists(file.path("..", filename))) {
  print(paste(filename, " has been created!!!", sep = ""))
} else{
  print(paste(filename, " cannot be created!!!", sep = ""))
}



# Print some empty space
cat(rep("\n", 3))
