#' Clean data for function sps()
#' @param X Site-level variables for the target population of sites. Row names should be names of sites. X cannot contain missing data.
#' @param scale TRUE or FALSE. Whether to standardize numeric variables to make each variable mean zero and standard deviation one.
#' @return \code{clean_for_sps} returns an object of \code{matrix} class, which we supply to \code{sps()}.
#'  \itemize{
#'    \item \code{X}: Site-level variables for the target population of sites.
#'  }
#' @references Egami and Lee. (2023+). Designing Multi-Context Studies for External Validity: Site Selection via Synthetic Purposive Sampling. Available at \url{https://naokiegami.com/paper/sps.pdf}.
#' @export
clean_for_sps <- function(X, scale = TRUE){

  ## Check
  if(any(is.na(X))){
    stop("'X' contains missing data. Please supply X without missing values. Consider using `impute_var()` in R package `spsRdata`. ")
  }

  ##
  if(is.null(rownames(X)) == TRUE){
    warning("'X' does not have row names. Please specify row names of 'X' ")
  }else if(all(rownames(X) == seq(1:nrow(X))) == TRUE){
    warning("rownames(X) is a sequence of numbers. Please consider specifying more informative row names for 'X' ")
  }

  ## Show all variables
  # cat(paste0(" ## 'X' includes the following variables: (", paste(colnames(X), collapse = ", "), ")"))

  ##
  class_var <- sapply(X, class)
  ind_non_num <- which(class_var == "character" | class_var == "factor")

  ## Check Variance of numeric variables
  if(any(class_var == "numeric")){
    var_X <- apply(X[, -ind_non_num, drop = FALSE], 2, var)
    if(any(var_X < 0.000001)){
      stop(" 'X' contains some variables with extremely small or no variation. Please check 'X'. For sps(), please only include variables to diversify in 'X'.  ")
    }
  }


  ## Transform non-numeric variables
  if(length(ind_non_num) > 0){
    unique_level <- apply(X[, ind_non_num, drop = FALSE], 2, function(x) length(unique(x)))
    if(any(unique_level > 0.9*nrow(X))){
      warning("Some non-numeric variables have many unique levels. For sps(), please only include variables to diversify. For example, 'X' should not include site names. ")
    }
  }

  if(length(ind_non_num) == 0){
    if(scale == TRUE){
      cat("## Standardizing numeric variables to make each variable mean zero and standard deviation one.\n")
      X_use <- base::scale(X)
    }
  }else if(length(ind_non_num) == 1){
    cat(paste0("## ", colnames(X)[ind_non_num], " is a non-numeric variable, and is binarized.\n"))
    X_non_num <- model.matrix(as.formula(paste0("~", colnames(X)[ind_non_num], " - 1")), data = as.data.frame(X))
    colnames(X_non_num) <- gsub(colnames(X)[ind_non_num], paste0(colnames(X)[ind_non_num], "_"), colnames(X_non_num))

    if(scale == TRUE){
      cat("## Standardizing numeric variables to make each variable mean zero and standard deviation one.\n")
      X_num <- base::scale(X[, -ind_non_num])
    }else{
      X_num <- X
    }
    X_use <- cbind(X_num, X_non_num)

  }else if(length(ind_non_num) > 1){
    cat(paste0("## (", paste(colnames(X)[ind_non_num], collapse = ", "), ") are non-numeric variables, and are binarized.\n"))
    X_non_num_list <- list()
    for(z in 1:length(ind_non_num)){
      X_non_num <- model.matrix(as.formula(paste0("~", colnames(X)[ind_non_num[z]], " - 1")), data = as.data.frame(X))
      colnames(X_non_num) <- gsub(colnames(X)[ind_non_num[z]], paste0(colnames(X)[ind_non_num[z]], "_"), colnames(X_non_num))
      X_non_num_list[[z]] <- X_non_num
    }
    X_non_num_use <- do.call("cbind", X_non_num_list)

    if(scale == TRUE){
      cat("## Standardizing numeric variables to make each variable mean zero and standard deviation one.\n")
      X_num <- base::scale(X[, -ind_non_num])
    }else{
      X_num <- X
    }
    X_use <- cbind(X_num, X_non_num_use)
  }

  ## check scale
  if(scale == FALSE){
    X_sd <- apply(X_use, 2, sd)
    if(max(X_sd)/min(X_sd) >= 100){
      warning("Some variables in X have standard deviation more than 100 times larger than other variables. This might cause estimation problem in sps(). Please consider using `scale()` to make standard deviations of variables comparable. ")
    }
  }

  return(X_use)

}
