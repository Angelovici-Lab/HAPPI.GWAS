
## Outlier removal function
outlier_removal <- function(dat = NULL, by_column = c(1, 2), start_column = 3){
  
  # Convert first column to character
  dat[,1] <- as.character(dat[,1])
  
  # Convert the rest columns to numeric
  for (i in 2:ncol(dat)) {
    dat[,i] <- as.numeric(dat[,i])
  }
  
  # Create lmer formula
  if (length(by_column) > 1) {
    termlabels <- c()
    for (i in 1:length(by_column)) {
      temp <- paste("(1|", colnames(dat)[i], ")", sep = "")
      termlabels <- c(termlabels, temp)
    }
  }
  
  # Calculate threshold
  threshold <- stats::qt(1-.05/(2*nrow(dat)), (nrow(dat)-3))
  
  # Find outlier
  outliers_residuals <- apply(dat[, start_column:ncol(dat)], 2, FUN = function(x){
    if (length(by_column) > 1) {
      lme <- lme4::lmer(formula = stats::reformulate(termlabels = termlabels, response = "x"), data = dat, REML=TRUE)
    } else if(length(by_column) == 1){
      lme <- stats::lm(x ~ 1, data = dat)
    }
    res <- stats::residuals(lme)
    H <- hatvalues(lme)
    sigma <- summary(lme)$sigm
    sres <- sapply(1:length(res), function(i) res[[i]]/(sigma*sqrt(1-H[[i]])))
    which(abs(sres) > threshold)
  })
  
  if(!identical(outliers_residuals, integer(0))){
    temp_outliers_residuals <- outliers_residuals
    
    outlier_dat <- data.frame()
    
    # Remove outliers
    for (i in 1:length(temp_outliers_residuals)) {
      rows <- match(temp_outliers_residuals[[i]], row.names(dat[names(temp_outliers_residuals)[i]]))
      columns <- which(grepl(names(temp_outliers_residuals)[i], colnames(dat)))
      
      temp_outliers_residuals[[i]] <- dat[rows, c(by_column, columns)]
      
      dat[rows, columns] <- NA
      
      if(nrow(temp_outliers_residuals[[i]]) > 0){
        if(nrow(outlier_dat) == 0){
          outlier_dat <- temp_outliers_residuals[[i]]
        }else{
          outlier_dat <- merge(outlier_dat, temp_outliers_residuals[[i]], 
                               by = intersect(colnames(outlier_dat), colnames(temp_outliers_residuals[[i]])), 
                               all=TRUE)
        }
      }
      
    }
    
    # Re-arrange first column
    for (i in 1:length(by_column)) {
      dat <- dat[order(as.numeric(gsub("[[:alpha:]]", "", dat[,i]))),]
    }
    
    # Re-arrange row names
    row.names(dat) <- seq(from = 1, to = nrow(dat), by = 1)
    
  }
    
  # Return data
  if(exists("dat") & exists("outlier_dat") & exists("outliers_residuals")){
    return(list("Outlier_removed_data" = dat, "Outlier_data" = outlier_dat, "Outliers_residuals" = outliers_residuals))
  } else{
    return(dat)
  }
  
}