#' Stratification for Synthetic Purposive Sampling
#' @param X Site-level variables for the target population of sites. Row names should be names of sites.
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

  col_use <- condition[[1]]

  if((condition[[2]] %in% c("larger than or equal to",
                          "smaller than or equal to",
                          "between")) == FALSE){
    stop("The second element of `condition` should be `larger than or equal to`, `smaller than or equal to`, or `between` ")
  }
  if((num_site[[1]] %in% c("at least", "at most")) == FALSE){
    stop("The first element of `num_site` should be `larger than or equal to` or `smaller than or equal to` ")
  }

  type_use <- case_when(condition[[2]] == "larger than or equal to" ~ "geq",
                        condition[[2]] == "smaller than or equal to" ~ "leq",
                        condition[[2]] == "between" ~ "between")
  # if(type_use == "between"){
  #   value <- c(condition[3], condition[4])
  # }else{
  #   value <- c(condition[3])
  # }
  value <- condition[[3]]

  # C_l <- list()
  # for(i in 1:length(columns)){
  #   col_use <- columns[i]
  #   C_l[[i]] <- stratify_base(X = X, columns = col_use, type = type_use, value = value)
  # }
  # C_use <- do.call("rbind", C_l)

  C_use <- stratify_base(X = X, columns = col_use, type = type_use, value = value)

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

