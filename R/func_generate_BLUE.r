#' generate BLUE function
#'
#' @description The goal of generate_BLUE is to run Best Linear Unbiased Estimations.
#' @importFrom foreach %dopar%
#' @importFrom magrittr %>%
#' @param dat An input dataset.
#' @param by_column The accession column.
#' @param start_column The start column index for traits.
#' @return blue or NULL if something missing.
#' @keywords BLUE, Pre-GWAS
#' @export
#'
generate_BLUE <- function(dat = NULL, by_column = c(1, 2), start_column = 3){

  # Convert first column to character
  dat[,1] <- as.character(dat[,1])

  # Convert the rest columns to numeric
  for (i in 2:ncol(dat)) {
    dat[,i] <- as.numeric(dat[,i])
  }



  #######################################################################
  ## Outlier Removal
  #######################################################################

  # Create lmer formula
  if (length(by_column) > 0) {
    termlabels <- colnames(dat)[1]
    for (i in 2:length(by_column)) {
      temp <- paste("(1|", colnames(dat)[i], ")", sep = "")
      termlabels <- c(termlabels, temp)
    }
  }

  # Calculate threshold
  threshold <- qt(1-.05/(2*nrow(dat)), (nrow(dat)-3))

  # Find outlier
  outliers_residuals <- list()
  for(i in start_column:ncol(dat)){
    lme <- lmer(formula = reformulate(termlabels = termlabels, response = colnames(dat)[i]), data = dat, REML=TRUE)

    res <- residuals(lme)
    H <- hatvalues(lme)
    sigma <- summary(lme)$sigm
    sres <- sapply(1:length(res), function(i) res[[i]]/(sigma*sqrt(1-H[[i]])))
    outliers_residuals[[colnames(dat)[i]]] <- which(abs(sres) > threshold)
  }

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

  outlier_removed_dat = dat

  #######################################################################
  ## Box-cox Transformation
  #######################################################################

  # Create lmer formula
  if (length(by_column) > 0) {
    termlabels <- colnames(dat)[1]
    for (i in 2:length(by_column)) {
      temp <- paste("(1|", colnames(dat)[i], ")", sep = "")
      termlabels <- c(termlabels, temp)
    }
  }

  names <- colnames(dat[,start_column:ncol(dat)])

  # run transformation for each trait
  transformed_out <- list()
  for(i in start_column:ncol(dat)){
    lme <- lmer(formula = reformulate(termlabels = termlabels, response = colnames(dat)[i]), data = dat, REML=TRUE)
    transformed_out[[colnames(dat)[i]]] <- car::powerTransform(lme, family="bcnPower", lambda=c(-2, 2))
  }

  lambda <- list()

  # put lambdas in a list
  for(i in names(transformed_out)) {
    # isolate the lambda for each column in dat saved in transformed_out
    lambda[[i]] <- (transformed_out[[i]]$lambda)
  }

  if(length(lambda) > 0){
    lambda <- as.data.frame(lambda)
    temp <- lambda

    for (i in 2:nrow(dat)) {
      temp[i,] <- temp[i-1,]
    }

    dat[,start_column:ncol(dat)] <- dat[,start_column:ncol(dat)]^temp
  }

  # Re-arrange first column
  for (i in 1:length(by_column)) {
    dat <- dat[order(as.numeric(gsub("[[:alpha:]]", "", dat[,i]))),]
  }

  # Re-arrange row names
  row.names(dat) <- seq(from = 1, to = nrow(dat), by = 1)

  boxcox_transformed_dat = dat

  #######################################################################
  ## generate BLUE
  #######################################################################

  # Create lmer formula
  if (length(by_column) > 0) {
    termlabels <- colnames(dat)[1]
    for (i in 2:length(by_column)) {
      temp <- paste("(1|", colnames(dat)[i], ")", sep = "")
      termlabels <- c(termlabels, temp)
    }
  }

  blue = data.frame(stringsAsFactors = FALSE)

  # fit the model
  for(i in start_column:ncol(dat)){

    dat[is.infinite(dat[,i]),i] = NA
    dat[is.nan(dat[,i]),i] = NA

    lme <- lmer(formula = reformulate(termlabels = termlabels, response = colnames(dat)[i]), data = dat, REML=TRUE)

    # estimate BLUP
    modelblue <- lme4::fixef(lme)
    modelblue[2:length(modelblue)] = modelblue[2:length(modelblue)] + summary(lme)$coefficients[1]
    modelblue = as.data.frame(modelblue, stringsAsFactors = FALSE)
    colnames(modelblue)[1] = colnames(dat)[i]

    modelblue = tibble::rownames_to_column(.data = modelblue, var = "Line")

    if(nrow(blue) == 0){
      blue = modelblue
    } else{
      blue = blue %>% dplyr::full_join(modelblue, by = "Line")
    }

  }

  blue[,1] = sub("Line", "", blue[,1])

  blue$Line[1] = dat$Line[which(!(unique(dat[,1]) %in% blue[,1]))]

  blue <- blue[order(as.numeric(gsub("[[:alpha:]]", "", blue[,1]))),]



  if(exists("blue")){
    return(
      list(
        "Outlier_removed_data" = outlier_removed_dat, "Outlier_data" = outlier_dat, "Outliers_residuals" = outliers_residuals,
        "Lambda_values" = lambda, "Boxcox_transformed_data" = boxcox_transformed_dat,
        "BLUE" = blue
      )
    )
  } else{
    return(-1)
  }

}


