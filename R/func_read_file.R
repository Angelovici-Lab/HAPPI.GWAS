#' Read .csv file or tab delimited text file
#'
#' @description The goal of read_file is to read .csv file or tab delimited text file.
#' @param file_path An csv file or tab delimited text file file path.
#' @param header TRUE if headers are in the file.
#' @return data frame.
#' @export
#'
read_file <- function(file_path = NULL, header = TRUE){
  if (!is.null(file_path)) {
    if (file.exists(file_path)) {
      if (endsWith(file_path, ".csv")) {
        dat <- tryCatch({
                          utils::read.csv(
                            file_path, header = header, stringsAsFactors = FALSE, check.names = FALSE
                          )
                        }, error = function(e) {
          return(NULL)
        })
      } else if (endsWith(file_path, ".txt")) {
        dat <- tryCatch({
                          utils::read.table(
                            file_path, header = header, sep = "\t", quote = "",
                            stringsAsFactors = FALSE, check.names = FALSE
                          )
                        }, error = function(e) {
          return(NULL)
        })
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