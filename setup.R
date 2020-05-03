#!/usr/bin/env Rscript

## Clean the workspace
rm(list = ls())

# Print some empty space
cat(rep("\n", 3))

# Create these folders
foldernames <- c("raw_data", "BLUP_BLUE", "reference_files", "output", "yaml")
for (i in 1:length(foldernames)) {
  temp <- file.path("..", foldernames[i])
  if(!dir.exists(temp)){
    dir.create(temp)
    if(dir.exists(temp)){
      print(paste0(foldernames[i], " folder has been created!!!"))
    } else{
      print(paste0(foldernames[i], " folder cannot be created!!!"))
    }
  } else{
    print(paste0(foldernames[i], " folder exists!!!"))
  }
}

# Create Arabidopsis360.yaml file
filename <- "Arabidopsis360.yaml"
dat <- readLines(con = filename)
writeLines(dat, con = file.path("..", filename))
if (file.exists(file.path("..", filename))) {
  print(paste0(filename, " has been created!!!"))
} else{
  print(paste0(filename, " cannot be created!!!"))
}

# Create Arabidopsis1001.yaml file 
filename <- "Arabidopsis1001.yaml"
dat <- readLines(con = filename)
writeLines(dat, con = file.path("..", filename))
if (file.exists(file.path("..", filename))) {
  print(paste0(filename, " has been created!!!"))
} else{
  print(paste0(filename, " cannot be created!!!"))
}

# Create Maize.yaml file
filename <- "Maize.yaml"
dat <- readLines(con = filename)
writeLines(dat, con = file.path("..", filename))
if (file.exists(file.path("..", filename))) {
  print(paste0(filename, " has been created!!!"))
} else{
  print(paste0(filename, " cannot be created!!!"))
}

# Print some empty space
cat(rep("\n", 3))
