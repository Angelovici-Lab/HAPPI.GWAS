
if(!require(yaml)) install.packages("yaml", dependencies = TRUE)
if(!require(argparse)) install.packages("argparse", dependencies = TRUE)
if(!require(doParallel)) install.packages("doParallel", dependencies = TRUE)

library(HAPPI.GWAS)

cores = 3
fp = file.path("/home/ycth8/data/projects/HAPPI_GWAS/HAPPI_GWAS/Demo_GLM.yaml")


## Read .csv file or tab delimited text file
read_file <- function(file_path = NULL, header = TRUE){
  if (!is.null(file_path)) {
    if (file.exists(file_path)) {
      if (endsWith(file_path, ".csv")) {
        dat <- try(utils::read.csv(file_path, header = header, stringsAsFactors = FALSE, check.names = FALSE))
      } else if (endsWith(file_path, ".txt")) {
        dat <- try(utils::read.table(file_path, header = header, sep = "\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE))
      } else{
        # unknown file extension.
        return(NULL)
      }
      return(dat)
    } else{
      # file_path does not exist.
    }
  } else{
    # file_path path is NULL.
  }

  return(NULL)
}



#######################################################################
## Starting and initialization
#######################################################################

cat(rep("\n", 2));print("-------------------- HAPPI_GWAS Start --------------------");cat(rep("\n", 2))

registerDoParallel(cores = cores)

cat(rep("\n", 2))
print(paste("Number of cores will be used is ", getDoParWorkers(), sep = ""))
cat(rep("\n", 2))

#######################################################################
## Read in YAML file data
#######################################################################

## Import configuration file and read in raw data and reference files
yaml_dat <- tryCatch({
  read_yaml(fp)
}, error = function(e) {
  print("The yaml file is invalid!!!")
  cat(rep("\n", 2))
  # print_help()
  quit(status = -1)
})

#######################################################################
## Check yaml data parameters
#######################################################################

if (exists("yaml_dat")) {
  if(is.null(yaml_dat$raw_data) & is.null(yaml_dat$BLUP)){
    print("The yaml data does not have input file path!!!")
    quit(status = -1)
  }
  if(is.null(yaml_dat$by_column) | is.null(yaml_dat$start_column) |
      is.null(yaml_dat$BLUP_by_column) | is.null(yaml_dat$BLUP_start_column)){
        print("The yaml data does not have input file by_column or start_column configuration!!!")
        quit(status = -1)
  }
  if(is.null(yaml_dat$output)){
    print("The yaml data does not have output path!!!")
    quit(status = -1)
  }
} else{
  print("The yaml data was not read into the work space!!!")
  quit(status = -1)
}



#######################################################################
## Read in all file and data that are specified in the YAML file
#######################################################################

cat(rep("\n", 2))
## Import raw data
raw_data <- read_file(file_path = yaml_dat$raw_data)
if (is.null(raw_data)) {
  print("The raw_data parameter is NULL.")
} else{
  print("raw_data has been loaded into memory.")
}

if (!is.null(yaml_dat$by_column)) {
  by_column <- yaml_dat$by_column
  print(paste("by_column: ", by_column, sep = ""))
} else{
  print("The by_column parameter is NULL.")
}

if (!is.null(yaml_dat$start_column)) {
  start_column <- yaml_dat$start_column
  print(paste("start_column: ", start_column, sep = ""))
} else{
  print("The start_column parameter is NULL.")
}

cat(rep("\n", 2))
## Import BULP data
BLUP <- read_file(file_path = yaml_dat$BLUP)
if (is.null(BLUP)) {
  print("The BLUP parameter is NULL.")
} else{
  print("BLUP has been loaded into memory.")
}

if (!is.null(yaml_dat$BLUP_by_column)) {
  BLUP_by_column <- yaml_dat$BLUP_by_column
  print(paste("BLUP_by_column: ", BLUP_by_column, sep = ""))
} else{
  print("The BLUP_by_column parameter is NULL.")
}

if (!is.null(yaml_dat$BLUP_start_column)) {
  BLUP_start_column <- yaml_dat$BLUP_start_column
  print(paste("BLUP_start_column: ", BLUP_start_column, sep = ""))
} else{
  print("The BLUP_start_column parameter is NULL.")
}

cat(rep("\n", 2))
## Import GAPIT reference files
# GAPIT kinship matrix
GAPIT_kinship_matrix <- read_file(file_path = yaml_dat$GAPIT_kinship_matrix, header = FALSE)
if (is.null(GAPIT_kinship_matrix)) {
  print("The GAPIT_kinship_matrix parameter is NULL.")
} else{
  print("GAPIT_kinship_matrix has been loaded into memory.")
}

# GAPIT covariates
GAPIT_covariates <- read_file(file_path = yaml_dat$GAPIT_covariates)
if (is.null(GAPIT_covariates)) {
  print("The GAPIT_covariates parameter is NULL.")
} else{
  print("GAPIT_covariates has been loaded into memory.")
}

# GAPIT hapmap file
GAPIT_hapmap <- read_file(file_path = yaml_dat$GAPIT_hapmap, header = FALSE)
if (is.null(GAPIT_hapmap)) {
  print("The GAPIT_hapmap parameter is NULL.")
} else{
  print("GAPIT_hapmap has been loaded into memory.")
}

# GAPIT genotype data (numeric)
GAPIT_genotype_data_numeric <- read_file(file_path = yaml_dat$GAPIT_genotype_data_numeric)
if (is.null(GAPIT_genotype_data_numeric)) {
  print("The GAPIT_genotype_data_numeric parameter is NULL.")
} else{
  print("GAPIT_genotype_data_numeric has been loaded into memory.")
}

# GAPIT genotype map (numeric)
GAPIT_genotype_map_numeric <- read_file(file_path = yaml_dat$GAPIT_genotype_map_numeric)
if (is.null(GAPIT_genotype_map_numeric)) {
  print("The GAPIT_genotype_map_numeric parameter is NULL.")
} else{
  print("GAPIT_genotype_map_numeric has been loaded into memory.")
}

# GAPIT hapmap file extension
if (!is.null(yaml_dat$GAPIT_hapmap_file_extension)) {
  GAPIT_hapmap_file_extension <- yaml_dat$GAPIT_hapmap_file_extension
  print(paste("GAPIT_hapmap_file_extension: ", GAPIT_hapmap_file_extension, sep = ""))
} else{
  GAPIT_hapmap_file_extension <- NULL
  print("The GAPIT_hapmap_file_extension parameter is NULL.")
}

# GAPIT genotype data numeric file extension
if (!is.null(yaml_dat$GAPIT_genotype_data_numeric_file_extension)) {
  GAPIT_genotype_data_numeric_file_extension <- yaml_dat$GAPIT_genotype_data_numeric_file_extension
  print(paste("GAPIT_genotype_data_numeric_file_extension: ", GAPIT_genotype_data_numeric_file_extension, sep = ""))
} else{
  GAPIT_genotype_data_numeric_file_extension <- NULL
  print("The GAPIT_genotype_data_numeric_file_extension parameter is NULL.")
}

# GAPIT genotype map numeric file extension
if (!is.null(yaml_dat$GAPIT_genotype_map_numeric_file_extension)) {
  GAPIT_genotype_map_numeric_file_extension <- yaml_dat$GAPIT_genotype_map_numeric_file_extension
  print(paste("GAPIT_genotype_map_numeric_file_extension: ", GAPIT_genotype_map_numeric_file_extension, sep = ""))
} else{
  GAPIT_genotype_map_numeric_file_extension <- NULL
  print("The GAPIT_genotype_map_numeric_file_extension parameter is NULL.")
}

# GAPIT hapmap filename
if (!is.null(yaml_dat$GAPIT_hapmap_filename)) {
  GAPIT_hapmap_filename <- yaml_dat$GAPIT_hapmap_filename
  print(paste("GAPIT_hapmap_filename: ", GAPIT_hapmap_filename, sep = ""))
} else{
  GAPIT_hapmap_filename <- NULL
  print("The GAPIT_hapmap_filename parameter is NULL.")
}

# GAPIT genotype data numeric filename
if (!is.null(yaml_dat$GAPIT_genotype_data_numeric_filename)) {
  GAPIT_genotype_data_numeric_filename <- yaml_dat$GAPIT_genotype_data_numeric_filename
  print(paste("GAPIT_genotype_data_numeric_filename: ", GAPIT_genotype_data_numeric_filename, sep = ""))
} else{
  GAPIT_genotype_data_numeric_filename <- NULL
  print("The GAPIT_genotype_data_numeric_filename parameter is NULL.")
}

# GAPIT genotype map numeric filename
if (!is.null(yaml_dat$GAPIT_genotype_map_numeric_filename)) {
  GAPIT_genotype_map_numeric_filename <- yaml_dat$GAPIT_genotype_map_numeric_filename
  print(paste("GAPIT_genotype_map_numeric_filename: ", GAPIT_genotype_map_numeric_filename, sep = ""))
} else{
  GAPIT_genotype_map_numeric_filename <- NULL
  print("The GAPIT_genotype_map_numeric_filename parameter is NULL.")
}

# GAPIT genotype file path
if (!is.null(yaml_dat$GAPIT_genotype_file_path)) {
  GAPIT_genotype_file_path <- file.path(yaml_dat$GAPIT_genotype_file_path)
  if(dir.exists(file.path(GAPIT_genotype_file_path))){
    print(paste("GAPIT_genotype_file_path: ", GAPIT_genotype_file_path, sep = ""))
  } else{
    print("The GAPIT_genotype_file_path parameter is a file path that does not exists.")
    quit(status = -1)
  }
} else{
  GAPIT_genotype_file_path <- NULL
  print("The GAPIT_genotype_file_path parameter is NULL.")
}

# GAPIT genotype file named sequentially from
if (!is.null(yaml_dat$GAPIT_genotype_file_named_sequentially_from)) {
  GAPIT_genotype_file_named_sequentially_from <- as.numeric(yaml_dat$GAPIT_genotype_file_named_sequentially_from)
  print(paste("GAPIT_genotype_file_named_sequentially_from: ", GAPIT_genotype_file_named_sequentially_from, sep = ""))
} else{
  GAPIT_genotype_file_named_sequentially_from <- 0
  print("The GAPIT_genotype_file_named_sequentially_from parameter is 0.")
}

# GAPIT genotype file named sequentially to
if (!is.null(yaml_dat$GAPIT_genotype_file_named_sequentially_to)) {
  GAPIT_genotype_file_named_sequentially_to <- as.numeric(yaml_dat$GAPIT_genotype_file_named_sequentially_to)
  print(paste("GAPIT_genotype_file_named_sequentially_to: ", GAPIT_genotype_file_named_sequentially_to, sep = ""))
} else{
  GAPIT_genotype_file_named_sequentially_to <- 0
  print("The GAPIT_genotype_file_named_sequentially_to parameter is 0.")
}

# GAPIT model
if (!is.null(yaml_dat$GAPIT_model)) {
  GAPIT_model <- yaml_dat$GAPIT_model
  print(paste("GAPIT_model: ", GAPIT_model, sep = ""))
} else{
  GAPIT_model <- NULL
  print("The GAPIT_model parameter is NULL.")
}

# GAPIT SNP.MAF
if (!is.null(yaml_dat$GAPIT_SNP_MAF)) {
  GAPIT_SNP_MAF <- as.numeric(yaml_dat$GAPIT_SNP_MAF)
  if(is.na(GAPIT_SNP_MAF)) { GAPIT_SNP_MAF <- 0 }
  print(paste("GAPIT_SNP_MAF: ", GAPIT_SNP_MAF, sep = ""))
} else{
  GAPIT_SNP_MAF <- 0
  print("The GAPIT_SNP_MAF parameter is 0.")
}

# GAPIT PCA total
if (!is.null(yaml_dat$GAPIT_PCA_total)) {
  GAPIT_PCA_total <- as.numeric(yaml_dat$GAPIT_PCA_total)
  if(is.na(GAPIT_PCA_total)) { GAPIT_PCA_total <- 0 }
  print(paste("GAPIT_PCA_total: ", GAPIT_PCA_total, sep = ""))
} else{
  GAPIT_PCA_total <- 0
  print("The GAPIT_PCA_total parameter is 0.")
}

# GAPIT Model selection
if (!is.null(yaml_dat$GAPIT_Model_selection)) {
  GAPIT_Model_selection <- as.logical(yaml_dat$GAPIT_Model_selection)
  if(!is.logical(GAPIT_Model_selection) | is.na(GAPIT_Model_selection)) { GAPIT_Model_selection <- FALSE }
  print(paste("GAPIT_Model_selection: ", GAPIT_Model_selection, sep = ""))
} else{
  GAPIT_Model_selection <- FALSE
  print("The GAPIT_Model_selection parameter is FALSE.")
}

# GAPIT SNP test
if (!is.null(yaml_dat$GAPIT_SNP_test)) {
  GAPIT_SNP_test <- as.logical(yaml_dat$GAPIT_SNP_test)
  if(!is.logical(GAPIT_SNP_test) | is.na(GAPIT_SNP_test)) { GAPIT_SNP_test <- FALSE }
  print(paste("GAPIT_SNP_test: ", GAPIT_SNP_test, sep = ""))
} else{
  GAPIT_SNP_test <- FALSE
  print("The GAPIT_SNP_test parameter is FALSE.")
}

# GAPIT file output
if (!is.null(yaml_dat$GAPIT_file_output)) {
  GAPIT_file_output <- as.logical(yaml_dat$GAPIT_file_output)
  if(!is.logical(GAPIT_file_output) | is.na(GAPIT_file_output)) { GAPIT_file_output <- FALSE }
  print(paste("GAPIT_file_output: ", GAPIT_file_output, sep = ""))
} else{
  GAPIT_file_output <- FALSE
  print("The GAPIT_file_output parameter is FALSE.")
}

# GAPIT p.Value threshold
if (!is.null(yaml_dat$GAPIT_p_value_threshold)) {
  GAPIT_p_value_threshold <- as.numeric(yaml_dat$GAPIT_p_value_threshold)
  if(is.logical(GAPIT_p_value_threshold)){ GAPIT_p_value_threshold <- NA }
  print(paste("GAPIT_p_value_threshold: ", GAPIT_p_value_threshold, sep = ""))
} else{
  GAPIT_p_value_threshold <- NA
  print("The GAPIT_p_value_threshold parameter is NA.")
}

# GAPIT p.Value FDR threshold
if (!is.null(yaml_dat$GAPIT_p_value_fdr_threshold)) {
  GAPIT_p_value_fdr_threshold <- as.numeric(yaml_dat$GAPIT_p_value_fdr_threshold)
  if(is.logical(GAPIT_p_value_fdr_threshold)){ GAPIT_p_value_fdr_threshold <- NA }
  print(paste("GAPIT_p_value_fdr_threshold: ", GAPIT_p_value_fdr_threshold, sep = ""))
} else{
  GAPIT_p_value_fdr_threshold <- NA
  print("The GAPIT_p_value_fdr_threshold parameter is NA.")
}

# GAPIT LD_number
if (!is.null(yaml_dat$GAPIT_LD_number)) {
  GAPIT_LD_number <- as.numeric(yaml_dat$GAPIT_LD_number)
  print(paste("GAPIT_LD_number: ", GAPIT_LD_number, sep = ""))
} else{
  print("The GAPIT_LD_number parameter is NULL.")
}

cat(rep("\n", 2))
## Haploview
# Haploview_file_path
if (!is.null(yaml_dat$Haploview_file_path)) {
  Haploview_file_path <- normalizePath(file.path(yaml_dat$Haploview_file_path))
  if(dir.exists(file.path(Haploview_file_path))){
    print(paste("Haploview_file_path: ", Haploview_file_path, sep = ""))
  } else{
    print("The Haploview_file_path parameter is a file path that does not exists.")
    quit(status = -1)
  }
} else{
  Haploview_file_path <- NULL
  print("The Haploview_file_path parameter is NULL.")
}

# Haploview_file_name
if (!is.null(yaml_dat$Haploview_file_name)) {
  Haploview_file_name <- as.character(yaml_dat$Haploview_file_name)
  print(paste("Haploview_file_name: ", Haploview_file_name, sep = ""))
} else{
  Haploview_file_name <- NULL
  print("The Haploview_file_name parameter is NULL.")
}

# Haploview_file_extension
if (!is.null(yaml_dat$Haploview_file_extension)) {
  Haploview_file_extension <- as.character(yaml_dat$Haploview_file_extension)
  print(paste("Haploview_file_extension: ", Haploview_file_extension, sep = ""))
} else{
  Haploview_file_extension <- NULL
  print("The Haploview_file_extension parameter is NULL.")
}

# Haploview_file_named_sequentially_from
if (!is.null(yaml_dat$Haploview_file_named_sequentially_from)) {
  Haploview_file_named_sequentially_from <- as.numeric(yaml_dat$Haploview_file_named_sequentially_from)
  print(paste("Haploview_file_named_sequentially_from: ", Haploview_file_named_sequentially_from, sep = ""))
} else{
  Haploview_file_named_sequentially_from <- NA
  print("The Haploview_file_named_sequentially_from parameter is NA.")
}

# Haploview_file_named_sequentially_to
if (!is.null(yaml_dat$Haploview_file_named_sequentially_to)) {
  Haploview_file_named_sequentially_to <- as.numeric(yaml_dat$Haploview_file_named_sequentially_to)
  print(paste("Haploview_file_named_sequentially_to: ", Haploview_file_named_sequentially_to, sep = ""))
} else{
  Haploview_file_named_sequentially_to <- NA
  print("The Haploview_file_named_sequentially_to parameter is NA.")
}

cat(rep("\n", 2))
## Match Gene Start and Gene Stop
# GFF_file_path
if (!is.null(yaml_dat$GFF_file_path)) {
  GFF_file_path <- normalizePath(file.path(yaml_dat$GFF_file_path))
  if(dir.exists(file.path(GFF_file_path))){
    print(paste("GFF_file_path: ", GFF_file_path, sep = ""))
  } else{
    print("The GFF_file_path parameter is a file path that does not exists.")
    quit(status = -1)
  }
} else{
  GFF_file_path <- NULL
  print("The GFF_file_path parameter is NULL.")
}

# GFF_file_name
if (!is.null(yaml_dat$GFF_file_name)) {
  GFF_file_name <- as.character(yaml_dat$GFF_file_name)
  print(paste("GFF_file_name: ", GFF_file_name, sep = ""))
} else{
  GFF_file_name <- NULL
  print("The GFF_file_name parameter is NULL.")
}

# GFF_file_extension
if (!is.null(yaml_dat$GFF_file_extension)) {
  GFF_file_extension <- as.character(yaml_dat$GFF_file_extension)
  print(paste("GFF_file_extension: ", GFF_file_extension, sep = ""))
} else{
  GFF_file_extension <- NULL
  print("The GFF_file_extension parameter is NULL.")
}

# GFF_file_named_sequentially_from
if (!is.null(yaml_dat$GFF_file_named_sequentially_from)) {
  GFF_file_named_sequentially_from <- as.numeric(yaml_dat$GFF_file_named_sequentially_from)
  print(paste("GFF_file_named_sequentially_from: ", GFF_file_named_sequentially_from, sep = ""))
} else{
  GFF_file_named_sequentially_from <- NA
  print("The GFF_file_named_sequentially_from parameter is NA.")
}

# GFF_file_named_sequentially_to
if (!is.null(yaml_dat$GFF_file_named_sequentially_to)) {
  GFF_file_named_sequentially_to <- as.numeric(yaml_dat$GFF_file_named_sequentially_to)
  print(paste("GFF_file_named_sequentially_to: ", GFF_file_named_sequentially_to, sep = ""))
} else{
  GFF_file_named_sequentially_to <- NA
  print("The GFF_file_named_sequentially_to parameter is NA.")
}

cat(rep("\n", 2))
## Create output folder
if (!is.null(yaml_dat$output)) {
  if (!dir.exists(file.path(yaml_dat$output))) {
    dir.create(path = file.path(yaml_dat$output), showWarnings = TRUE, recursive = TRUE)
    if (dir.exists(file.path(yaml_dat$output))) {
      output <- normalizePath(file.path(yaml_dat$output))
      print(paste0("The output folder is created. Output path: ", output))
    } else{
      print("The output folder cannot be created.")
      quit(status = -1)
    }
  } else{
    output <- normalizePath(file.path(yaml_dat$output))
    print(paste0("The output folder exists. Output path: ", output))
  }
} else{
  print("The output parameter is NULL.")
  quit(status = -1)
}

cat(rep("\n", 2))

#######################################################################
## Run action base on arguments
#######################################################################

if (exists("BLUP") & !is.null(BLUP) & exists("BLUP_by_column") & exists("BLUP_start_column") &
  exists("GAPIT_LD_number") & GAPIT_LD_number >= 0 & dir.exists(output)) {

    folder_path <- file.path(output, "GAPIT")

    if (!dir.exists(folder_path)) {
      dir.create( path = folder_path, showWarnings = TRUE, recursive = TRUE)
    } else{
      print("The GAPIT folder exists.")
    }

    # Using customized function to run GAPIT
    combined_gwas_result <-
      farming_with_GAPIT(
        dat = BLUP,
        by_column = BLUP_by_column,
        start_column = BLUP_start_column,
        output_path = folder_path,
        p_value_threshold = GAPIT_p_value_threshold,
        p_value_fdr_threshold = GAPIT_p_value_fdr_threshold,
        ld_number = GAPIT_LD_number,
        KI = GAPIT_kinship_matrix,
        CV = GAPIT_covariates,
        G = GAPIT_hapmap,
        GD = GAPIT_genotype_data_numeric,
        GM = GAPIT_genotype_map_numeric,
        file.Ext.G = GAPIT_hapmap_file_extension,
        file.Ext.GD = GAPIT_genotype_data_numeric_file_extension,
        file.Ext.GM = GAPIT_genotype_map_numeric_file_extension,
        file.G = GAPIT_hapmap_filename,
        file.GD = GAPIT_genotype_data_numeric_filename,
        file.GM = GAPIT_genotype_map_numeric_filename,
        file.path = GAPIT_genotype_file_path,
        file.from = GAPIT_genotype_file_named_sequentially_from,
        file.to = GAPIT_genotype_file_named_sequentially_to,
        model = GAPIT_model,
        SNP.MAF = GAPIT_SNP_MAF,
        PCA.total = GAPIT_PCA_total,
        Model.selection = GAPIT_Model_selection,
        SNP.test = GAPIT_SNP_test,
        file.output = GAPIT_file_output
      )

}

# if combined_gwas_result is not null and other requirements are satisfied then extract haplotype
if(exists("combined_gwas_result") & !is.null(combined_gwas_result) & !is.null(Haploview_file_path) & !is.null(Haploview_file_name) &
  !is.null(Haploview_file_extension) & !is.na(Haploview_file_named_sequentially_from) & !is.na(Haploview_file_named_sequentially_to)){

    folder_path <- file.path(output, "GAPIT")

    if (!dir.exists(folder_path)) {
      dir.create( path = folder_path, showWarnings = TRUE, recursive = TRUE)
    } else{
      print("The GAPIT folder exists.")
    }

    combined_gwas_result <- extract_haplotype(
      combined_gwas_result = combined_gwas_result,
      output_path = folder_path,
      Haploview_file_path = Haploview_file_path,
      Haploview_file_name = Haploview_file_name,
      Haploview_file_extension = Haploview_file_extension,
      Haploview_file_named_sequentially_from = Haploview_file_named_sequentially_from,
      Haploview_file_named_sequentially_to = Haploview_file_named_sequentially_to
    )
}

# if combined_gwas_result is not null and other requirements are satisfied then search genes
if(exists("combined_gwas_result") & !is.null(combined_gwas_result) & !is.null(GFF_file_path) & !is.null(GFF_file_name) &
  !is.null(GFF_file_extension) & !is.na(GFF_file_named_sequentially_from) & !is.na(GFF_file_named_sequentially_to)){

    folder_path <- file.path(output, "GAPIT")

    if (!dir.exists(folder_path)) {
      dir.create( path = folder_path, showWarnings = TRUE, recursive = TRUE)
    } else{
      print("The GAPIT folder exists.")
    }

    combined_gwas_result <- search_genes(
      combined_gwas_result = combined_gwas_result,
      output_path = folder_path,
      GFF_file_path = GFF_file_path,
      GFF_file_name = GFF_file_name,
      GFF_file_extension = GFF_file_extension,
      GFF_file_named_sequentially_from = GFF_file_named_sequentially_from,
      GFF_file_named_sequentially_to = GFF_file_named_sequentially_to
    )
}
