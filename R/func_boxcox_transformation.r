
## Box cox transformation function
boxcox_transformation <- function(dat = NULL, by_column = c(1, 2), start_column = 3){
    if(!is.null(dat)){
        # Convert first column to character
        dat[,1] <- as.character(dat[,1])

        # Convert the rest columns to numeric
        for (i in 2:ncol(dat)) {
            dat[,i] <- as.numeric(dat[,i])
        }

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
        transformed_out <- apply(dat[,start_column:ncol(dat)], 2, FUN = function(x){
            lme <- lmer(formula = reformulate(termlabels = termlabels, response = "x"), data = dat, REML=TRUE)
            car::powerTransform(lme, family="bcPower", lambda=c(-2, 2))
        })

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

        if(exists("lambda") & length(lambda) > 0 & exists("dat")) {
            return(list("Lambda_values" = lambda, "Boxcox_transformed_data" = dat))
        } else{
            return(dat)
        }
    } else{
        return(NULL)
    }
  
}
