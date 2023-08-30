#' Stratification for Synthetic Purposive Sampling
#' @param X Site-level variables for the target population of sites. Row names should be names of sites. X cannot contain missing data.
#' @param num_site A list of two elements, e.g., \code{list("at least", 1)}. This argument specifies the number of sites that should satisfy \code{condition} specified below. The first element should be either \code{at least} or \code{at most}. The second element is integer. For example, \code{list("at least", 1)} means that we stratify SPS such that we select *at least 1* site that satisfies \code{condition} (specified below).
#' @param condition A list of three elements, e.g., \code{list("GDP", "larger than or equal to", 1)}. This argument specifies conditions for stratification. The first element should be a name of a site-level variable. The second element should be either \code{larger than or equal to}, \code{smaller than or equal to}, or \code{between}. The third element is a vector of length 1 or 2. When the second element is \code{between}, the third element should be a vector of two values. For example, \code{list("GDP", "larger than or equal to", 1)} means that we stratify SPS such that we select \code{num_site} sites that have *GDP larger than or equal to 1*.
#' @return \code{stratify_sps} returns an object of \code{stratify_sps} class, which we supply to \code{sps()}.
#'  \itemize{
#'    \item \code{C}: A matrix on the left-hand side of linear constraints. The number of columns is the number of sites in the target population (=\code{nrow(X)}) and the number of rows is the number of constraints.
#'    \item \code{c0}: A vector on the right-hand side of linear constraints. The length is the number of constraints.
#'  }
#' @references Egami and Lee. (2023+). Designing Multi-Context Studies for External Validity: Site Selection via Synthetic Purposive Sampling. Available at \url{https://naokiegami.com/paper/sps.pdf}.
#' @export
stratify_sps <- function(X, num_site = NULL, condition = NULL){

  # ###############
  # Housekeeping
  # ###############
  ## X
  ## X contains NA
  if(any(is.na(X))){
    stop(" `X` contains missing data. Please supply X without missing values. Consider using `impute_var()` in R package `spsRdata`. ")
  }

  ## factor or character
  if(any(sapply(X, class) %in% c("factor", "character")) == TRUE){
    stop(" `X` contains `factor` or `character` variables. Before using sps(), please convert them into numeric or binary variables. ")
  }

  ## num_site
  if(is.list(num_site) == FALSE){
    stop(" `num_site` should be a list of two elements. ")
  }
  if((num_site[[1]] %in% c("at least", "at most")) == FALSE){
    stop(" The first element of `num_site` should be either `at least` or `at most` ")
  }
  if(is.numeric(num_site[[2]]) == FALSE){
    stop(" The second element of `num_site` should be class `numeric` ")
  }
  ## condition
  if(is.list(condition) == FALSE){
    stop(" `condition` should be a list of three elements. ")
  }
  if((condition[[1]] %in% colnames(X)) == FALSE){
    stop(" The first element of `condition` should be one of `colnames(X)` ")
  }
  if((condition[[2]] %in% c("larger than or equal to",
                            "smaller than or equal to",
                            "between")) == FALSE){
    stop(" The second element of `condition` should be `larger than or equal to`, `smaller than or equal to`, or `between` ")
  }
  if(is.numeric(condition[[3]]) == FALSE){
    stop(" The third element of `condition` should be class `numeric` ")
  }
  if(condition[[2]] == "between"){
    if(length(condition[[3]]) !=2){
      stop(" When the second element of `condition` is `between`, the third element of `condition` should be a vector of two values. ")
    }
  }
  ######

  # Creating C (left-hand-side)
  col_use <- condition[[1]]
  type_use <- case_when(condition[[2]] == "larger than or equal to" ~ "geq",
                        condition[[2]] == "smaller than or equal to" ~ "leq",
                        condition[[2]] == "between" ~ "between")
  value <- condition[[3]]
  C_use <- stratify_base(X = X, columns = col_use, type = type_use, value = value)

  if(num_site[[2]] %in% c(0, 1)){
    num_site_after <- "site"
  }else{
    num_site_after <- "sites"
  }
  if(sum(C_use) %in% c(0, 1)){
    sum_site_after <- "site satisfies"
  }else{
    sum_site_after <- "sites satisfy"
  }

  cat(paste0(sum(C_use), " ", sum_site_after, " the specified `condition` and sps() will select ",
             num_site[[1]], " ", num_site[[2]], " ", num_site_after, " from them.\n"))

  if(sum(C_use) == 0){
    warnings(paste0("There is no site that satisfies the specified `condition`"))
  }

  # Creating c0 (right-hand-side)
  if(num_site[[1]] == "at least"){
    c0_value <- as.numeric(num_site[[2]])
  }else if(num_site[1] == "at most"){
    c0_value <- -1*as.numeric(num_site[[2]])
  }
  c0_use <- rep(c0_value, nrow(C_use))

  st <- list("C" = C_use, "c0" = c0_use)
  class(st) <- c(class(st), "stratify_sps")
  return(st)
}

stratify_base <- function(X, columns = NULL, type = "geq", value = 1){

  if(is.null(columns)){
    columns <- seq(1:ncol(X))
  }
  C_mat <- matrix(NA, nrow = length(columns), ncol = nrow(X))
  for(i in 1:length(columns)){
    if(type == "geq"){
      C_mat[i,1:ncol(C_mat)] <- as.numeric(X[, columns[i]] >= value)
    }else if(type == "leq"){
      C_mat[i,1:ncol(C_mat)] <- as.numeric(X[, columns[i]] <= value)
    }else if(type == "between"){
      C_mat[i, 1:ncol(C_mat)] <- as.numeric(X[, columns[i]] >= value[1])*as.numeric(X[, columns[i]] <= value[2])
    }
  }
  return(C_mat)
}

