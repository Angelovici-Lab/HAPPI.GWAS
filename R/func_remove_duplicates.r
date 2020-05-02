
remove_duplicates <- function(dat = NULL, by_column = c(1,2)){
  
  temp <- list()
  
  for (i in 1:length(by_column)) {
    temp[[i]] <- dat[,i]
  }
  
  # Aggregate data, calculate means for duplicates
  dat <- stats::aggregate(x = dat, by = temp, FUN = "mean")
  
  for (i in 1:length(by_column)) {
    dat[,i+length(by_column)] <- dat[,i]
  }
  
  # Remove the aggregation grouping
  dat <- dat[, (by_column*-1)]
  
  # Re-arrange first column
  for (i in 1:length(by_column)) {
    dat <- dat[order(as.numeric(gsub("[[:alpha:]]", "", dat[,i]))),]
  }
  
  # Re-arrange row names
  row.names(dat) <- seq(from = 1, to = nrow(dat), by = 1)
  
  return(dat)
}
