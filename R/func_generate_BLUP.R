#' generate BLUP function
#'
#' @description The goal of generate_BLUP is to run Best Linear Unbiased Predictions.
#' @importFrom foreach %dopar%
#' @importFrom magrittr %>%
#' @importFrom data.table as.data.table
#' @param dat An input dataset.
#' @param by_column The accession column.
#' @param start_column The start column index for traits.
#' @return blue or NULL if something missing.
#' @keywords BLUP, Pre-GWAS
#' @export
#'
generate_BLUP <- function(dat = NULL, by_column = c(1, 2), start_column = 3){

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
    termlabels <- c()
    for (i in 1:length(by_column)) {
      temp <- paste("(1|", colnames(dat)[i], ")", sep = "")
      termlabels <- c(termlabels, temp)
    }
  }

  # Calculate threshold
  threshold <- qt(1-.05/(2*nrow(dat)), (nrow(dat)-3))

  # Find outlier
  outliers_residuals <- list()
  for(i in start_column:ncol(dat)){
    lme <- lme4::lmer(formula = reformulate(termlabels = termlabels, response = colnames(dat)[i]), data = dat, REML=TRUE)

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
    termlabels <- c()
    for (i in 1:length(by_column)) {
      temp <- paste("(1|", colnames(dat)[i], ")", sep = "")
      termlabels <- c(termlabels, temp)
    }
  }

  names <- colnames(dat[,start_column:ncol(dat)])

  # run transformation for each trait
  transformed_out <- list()
  for(i in start_column:ncol(dat)){
    lme <- lme4::lmer(formula = reformulate(termlabels = termlabels, response = colnames(dat)[i]), data = dat, REML=TRUE)
    transformed_out[[colnames(dat)[i]]] <- tryCatch({
                                                      car::powerTransform(lme, family="bcPower", lambda=c(-2, 2))
                                                    }, error = function (e){
      cat(rep("\n", 2))
      print(paste0("Lambda cannot be calculated for ", colnames(dat)[i]))
      return(1)
    })
  }

  lambda <- list()

  # put lambdas in a list
  for(i in names(transformed_out)) {
    # isolate the lambda for each column in dat saved in transformed_out
    lambda[[i]] <- tryCatch({
                              transformed_out[[i]]$lambda
                            },error = function (e) {
      return(transformed_out[[i]][1])
    })
  }

  if(length(lambda) > 0){
    lambda <- as.data.frame(lambda)
    lambda_matrix <- lambda

    for (i in 2:nrow(dat)) {
      lambda_matrix[i,] <- lambda_matrix[i-1,]
    }

    dat[,start_column:ncol(dat)] <- ifelse(
      lambda_matrix==0,
      log(dat[,start_column:ncol(dat)]),
      dat[,start_column:ncol(dat)]^lambda_matrix
    )
  }

  # Re-arrange first column
  for (i in 1:length(by_column)) {
    dat <- dat[order(as.numeric(gsub("[[:alpha:]]", "", dat[,i]))),]
  }

  # Re-arrange row names
  row.names(dat) <- seq(from = 1, to = nrow(dat), by = 1)

  not_transform_columns <- colnames(lambda)[lambda==1]
  lambda[lambda==1] <- NA

  boxcox_transformed_dat <- dat

  #######################################################################
  ## generate BLUP
  #######################################################################

  # Create lmer formula
  if (length(by_column) > 0) {
    termlabels <- c()
    for (i in 1:length(by_column)) {
      temp <- paste("(1|", colnames(dat)[i], ")", sep = "")
      termlabels <- c(termlabels, temp)
    }
  }

  # fit the model
  BLUP_out <- list()
  for(i in start_column:ncol(dat)){

    dat[is.infinite(dat[,i]),i] = NA
    dat[is.nan(dat[,i]),i] = NA

    lme <- lme4::lmer(formula = reformulate(termlabels = termlabels, response = colnames(dat)[i]), data = dat, REML=TRUE)

    # estimate BLUP
    modelblup <- lme4::ranef(lme)

    # extract BLUP and add grand mean when only one repetition present
    BLUP_out[[colnames(dat)[i]]] <- modelblup[[1]] + summary(lme)$coefficients[1]
  }

  if(length(BLUP_out) > 0){
    for (i in 1:length(BLUP_out)) {
      temp <- as.data.frame(BLUP_out[[i]])
      temp$New <- row.names(temp)
      temp$id <- names(BLUP_out)[i]
      temp <- temp[,c(3, 2, 1)]
      colnames(temp) <- c("id", "Line", "Intercept")
      row.names(temp) <- seq(from = 1, to = nrow(temp), by = 1)
      BLUP_out[[i]] <- temp

      if(i==1){
        BLUP_out_df <- BLUP_out[[i]]
      } else{
        BLUP_out_df <- rbind(BLUP_out_df, BLUP_out[[i]])
      }
    }

    blup <- reshape2::dcast(BLUP_out_df, Line ~ id, value.var="Intercept")

    colnames(blup)[1] <- colnames(dat)[1]
    blup <- blup[order(as.numeric(gsub("[[:alpha:]]", "", blup[,1]))),]
    row.names(blup) <- seq(from = 1, to = nrow(blup), by = 1)

  }

  not_converge_columns = c()
  for(i in 1:ncol(blup)){
    if(is.numeric(blup[,i])){
      if(var(blup[,i], na.rm = TRUE) < 1e-6){
        not_converge_columns = c(not_converge_columns, colnames(blup)[i])
      }
    }
  }
  if(!is.null(not_converge_columns)){
    blup = blup[, -(match(not_converge_columns, colnames(blup)))]
  }

  if(exists("blup")){
    return(
      list(
        "Outlier_removed_data" = outlier_removed_dat,
        "Outlier_data" = outlier_dat,
        "Outliers_residuals" = outliers_residuals,
        "Lambda_values" = lambda,
        "Boxcox_transformed_data" = boxcox_transformed_dat,
        "BLUP" = blup,
        "Not_converge_columns" = not_converge_columns,
        "Not_transform_columns" = not_transform_columns
      )
    )
  } else{
    return(-1)
  }

}
