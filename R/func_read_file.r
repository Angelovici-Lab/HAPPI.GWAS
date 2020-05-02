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
