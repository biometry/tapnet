#' Helper functions for tapnet
#'
#' Lower-level, non-exported functions to be called by the main tapnet functions
#'
#' They do roughly the following: 
#' \describe{
#'   \item{\code{bjornloglik}}{computes likelihood of obtaining the observed interaction matrix, given some expected matrix of interaction propabilities (P). In contrast to the multinomial, this function assumes that the marginal totals of P constrain the probabilities. Hence, if all observations for one column (or row) have been evaluated for their probabilities, any new observation cannot come from this column (or row) anymore. This means, the probabilities in that "depleted" column (or row) have to be proportionally distributed over the other rows (or columns). There are ongoing discussions among the authors, when this is the right approach. Function written by Björn Reineking, ISRAE Grenoble (many thanks!), hence the name.}
#'   \item{\code{fit_abund}}{computes expected P-matrix of observed interactions based only on abundances; if no external abundances are provided in the tapnet-object, it will use the marginal totals of the network(s) instead.}
#'   \item{\code{pems_from_tree}}{computes phylogenetic eigenvectors from a phylogenetic tree;}
#'   \item{\code{select_relevant_pems}}{identifies those phylogenetic eigenvectors (PEMs) of the full tree most relevant for a network containing only a subset of species;}
#'   \item{\code{tmatch}}{calculates interaction probabilities based on trait matching;}
#'   \item{\code{param_vec2list}}{converts a vector of parameters (for trait matching and latent trait combinations) into a named list;}
#'   \item{\code{loglik_tapnet}}{the log-likelihood function for fitting the tapnet model; actually quite an important function, easy to break, so not for the user to easily access;}
#'   \item{\code{latent_cor}}{computes correlation of fitted latent with true constructed traits for simulated data;}
#'   \item{\code{web_indices}}{computes the specified network indices for the provided network, after turning the prediction vector into a matrix;}
#'   \item{\code{refit_params}}{Simulate new networks from a fitted tapnet object, re-fit on the simulated network and output the parameter values.}
#'}
#'
#' @aliases helper functions
#' @aliases pems_from_tree select_relevant_pems tmatch param_vec2list loglik_tapnet latent_cor web_indices refit_params
#' 
#' @param tree phylogenetic tree in phylo format;
#' @param species a named vector of species, representing (some of) the tips of the phylogenetic tree;
#' @param delta_t vector of pairwise trait differences (higher - lower);
#' @param type trait matching function: either "normal" or "shiftlnorm";
#' @param width width parameter of trait matching function, similar to sd in the normal;
#' @param shift shift parameter (optimum trait distance), currently ignored in fitting;
#' @param err "baseline" probability of match, even if traits do not match at all; 
#' @param params  parameter vector with setting for tapnet simulation;
#' @param n  number of latent trait linear combination parameters (lower level);
#' @param m  number of latent trait linear combination parameters (higher level);
#' @param fit.delta logical; should the trait-weighting exponent delta be fitted?
#' @param networks the "networks" part of a tapnet object;
#' @param tmatch_type_pem type of trait matching function for latent traits;
#' @param tmatch_type_obs type(s) of trait matching functions for observed traits; can be a vector as long as there are traits;
#' @param TmatchMatrixList list of independent trait-matching matrices (one per network);
#' @param lambda LASSO shrinkage parameter to avoid collinearity issues when observed traits are phylogenetically correlated;
#' @param obj_function objective function, either "multinom" or "least squares" (leads to OLS fitting) or "bjorn" (leading to use of a somewhat weird but in some opinion the correct way to compute the likelihood);
#' @param true_pars parameters used for simulating the network;
#' @param fitted_pars parameters estimated by \code{\link{fit_tapnet}};
#' @param P matrix of same size as observed web, summing to 1, representing something like the probability of selecting a link; typically constructed as part of tapnet, i.e. the I-mat of fit_tapnet;
#' @param pems_low phylogenetic eigenvectors for the lower trophic level;
#' @param pems_high phylogenetic eigenvectors for the higher trophic level;
#' @param web data for an interaction network in vector form, e.g. from predict_tapnet;
#' @param web_dim vector of two numbers indicating the rows and column of the network matrix;
#' @param indices vector of names of network indices to compute; see \code{\link[bipartite]{networklevel}} for what is available;
#' @param tapnet a tapnet object;
#' @param fit a fitted tapnet;
#' @param fitted_I_mat the fitted I-matrix of a fitted tapnet object.
#' @param xi penalty for increasing the width of the (unstandardised) trait-matching function; defaults to 1 (implying a quadratically increasing weight of width in the optimisatino function, i.e. strongly acting against large values of sigma)

#'
#' @references Benadi et al. in prep
#'
#' @author Gita Benadi <gita.benadi@biom.uni-freiburg.de>, Carsten Dormann <carsten.dormann@biom.uni-freiburg.de> and Jochen Fründ <jochen.fruend@biom.uni-freiburg.de>
#'
#'
internalFunctions <- function(){}# just here to display the name for the help page through Roxygen correctly; grrrr

#' @rdname internalFunctions
bjornloglik <- function(web, P) {
  # Calculates likelihood of observed interactions in interaction matrix M
  # given expected probabilities P.
  # 
  # web: interaction matrix
  # P: probability matrix of the same dimensions as M
  
  if(any(dim(web) != dim(P))) stop("Dimensions of M and P do not match.")
  if(sum(P) != 1) {
    warning("P does not sum to 1. Will be normalised.")
    P <- P/sum(P, na.rm=TRUE)
  }
  loglik <- 0
  R <- rowSums(web)
  C <- colSums(web)
  for(i in 1:NROW(web)) {
    for(j in 1:NCOL(web)) {
      k <- web[i, j]
      if (k > 0) {
        loglik <- loglik + k * log(P[i, j])
        R[i] <- R[i] - k
        C[j] <- C[j] - k
        if(R[i] == 0) {
          r <- P[i,]
          c <- P[,j]
          P <- sweep(P, 2, 1-r, "/") # renormalise probabilities to account for "lost" interactions
          if(C[j] == 0) {
            P <- sweep(P, 1, 1-c, "/")
          }
        }
      }
    }
  }
  loglik
}

#' @rdname internalFunctions
fit_abund <- function(tapnet){
  # outputs P-matrix based only abundances
  Nwebs <- length(tapnet$networks)
  Ps <- list(1:Nwebs)
  for (i in 1:Nwebs){ # loop through all networks
    
    if (is.null(tapnet$networks[[i]]$abuns)){
      warning("Tapnet object contains no external abundances; will use matrix marginal totals instead.")
      theweb <- tapnet$networks[[i]]$web
      abundances <- list(rowSums(theweb), colSums(theweb))
    } else{
      abundances <- tapnet$networks[[i]]$abuns
    }
    # A <- tcrossprod(abundances[[1]], abundances[[2]]) # removes names
    A <- abundances[[1]] %*% t(abundances[[2]]) # partially keeps names
    rownames(A) <- names(abundances[[1]])
    Ps[[i]] <- A/sum(A)
  }
  return(Ps)
}
#fit_abund(tapnet_web1)


#' @rdname internalFunctions
pems_from_tree <- function(tree) {
  dirgraph <- Phylo2DirectedGraph(tree)
  pem <- as.data.frame(PEM.build(dirgraph))
  # sort alphabetically:
  pem <- pem[order(rownames(pem)),]
  return(pem)
}

#' @rdname internalFunctions
select_relevant_pems <- function(tree, # a phylogenetic tree
                                 species # species for which relevant PEMs should be extracted
                                 # (must be a subset of the tree's tip.label vector)
){
  # PEMs are computed from the full phylogentic information provided (using pems_from_tree);
  # for a given network, only some species may be present;
  # instead of computing PEMs only for those, we want to use the PEMs from the full tree, 
  # but only those that are useful for the species in the network;
  # this function computes the PEMs for the species present, then correlates them with the 
  # full set of PEMs and returns those PEMs that best correlate with the subtree's PEMs
  
  # check that all selected species appear in the tree:
  if (!all(species %in% tree$tip.label)) stop("Some selected species do not appear in the tree!")
  
  PEMsall <- pems_from_tree(tree)
  # keep only those species also present in the subweb:
  PEMsallsub <- PEMsall[rownames(PEMsall) %in% species, ]
  #reduce tree to those species present in the web:
  sub.tree <- keep.tip(tree, species)
  PEMssubtree <- pems_from_tree(sub.tree)
  # now correlate each subtreePEM with all PEMsall and identify the best correlation
  keep.these.PEMs <- NULL
  for (i in 1:ncol(PEMssubtree)){
    keep.these.PEMs[i] <- which.max(cor(cbind(PEMssubtree[,i], PEMsallsub))[-1,1])
  }
  # Note that this can be substantially fewer than the number of species
  # if several PEMs of the subweb are correlated with the same full web's PEM.
  out <- PEMsallsub[, sort(unique(keep.these.PEMs))]
  out <- out[order(rownames(out)),]
  return(out)
}

#' @rdname internalFunctions
tmatch <- function(delta_t, # Vector of pairwise trait differences (higher - lower)
                   type = "normal", # Trait matching function
                   width = 1, # Width parameter of trait matching function,
                   shift = 0, # shift parameter (optimum trait distance)
                   xi = 1, # weight that determines importance of penalty to value of width
                   err = 1E-5 # "baseline" probability of match, even if traits do not match at all
){# Calculate interaction probabilities based on trait matching
  # shift
  
  if (!(type %in% c("normal", "shiftlnormal"))) stop("Please use 'normal' or 'shiftlnorm' [in tmatch].")
  
  delta_t <- delta_t + shift
  
  # lognormal distribution with mode shifted to zero:
  if (type == "shiftlnorm") out <- dlnorm(delta_t + width, meanlog = log(width) + 1) + err
  
  # normal distribution, unstandardised but with sigma-based penalty:
  if (type == "normal") out <- 
    (sqrt(2*pi)*width) * dnorm(delta_t, mean = 0, sd = width) + err + xi * width^2 # penalty increases quadratically with sigma
    #(sqrt(2*pi)*width) * dnorm(delta_t, mean = 0, sd = width) + xi * dnorm(delta_t, mean = 0, sd = width) + err # penalty decreases linearly with sigma
  
  return(out)
}

#' @rdname internalFunctions
param_vec2list <- function(params, # Parameter vector
                           n, # Number of latent trait linear combination parameters (lower level)
                           m, # Number of latent trait linear combination parameters (higher level)
                           fit.delta=FALSE) {
  # Convert a vector of parameters (for trait matching and latent trait combinations)
  # to a named list
  par_list <- list()
  if (n > 0 & m > 0){ # first sort out the latents:
    par_list$lat_low <- params[1:n] # PEM linear combination parameters (lower level)
    par_list$lat_low[1] <- exp(params[1]) # Force positive value
    par_list$lat_high <- params[(n + 1):(n + m)] # PEM linear combination parameters (higher level)
    par_list$pem_shift <- params[n + m + 1] # PEM shift
    par_list$tmatch_width_pem <- exp(params[n + m + 2]) # PEM trait matching
  }
  
  if (fit.delta){ # if delta is to be fitted, it is the last of params
    if (length(params) > (n + m + 2)) {
      par_list$tmatch_width_obs <- exp(params[(n + m + 3):(length(params)-1)]) # observed traits matching
    } else { # no pems but observed traits:
      par_list$tmatch_width_obs <- exp(params[2:length(params)])
    }
    par_list$delta <- plogis(params[length(params)])
  } else { # if delta is NOT fitted, all last parameters go to traits
    if (length(params) > (n + m + 2)) {
      par_list$tmatch_width_obs <- exp(params[(n + m + 3):(length(params))]) # observed traits matching
    } else {
      par_list$tmatch_width_obs <- exp(params[2:length(params)])
    }
  }
  return(par_list)
}

#' @rdname internalFunctions
loglik_tapnet <- function(params, # Parameters (a *named* vector)
                          networks, # the "networks" part of a tapnet object
                          tmatch_type_pem, # Type of trait matching function for latent traits
                          tmatch_type_obs, # Type(s) of trait matching functions for observed traits
                          TmatchMatrixList=NULL, # independent trait-matching matrices (one per network)
                          lambda=0, # LASSO shrinkage parameter
                          obj_function = "multinom",  # Objective function: "least squares", "bjorn"
                          fit.delta=TRUE # either "multinom" (negative log-likelihood based on a multinomial distribution) or
                          # "sq_diff" (sum of squared differences between observed and predicted log(interactions + 1))
) {
  # Compute value of the objective function for fitting given the parameter values and "tap" data
  
  if (is.null(names(params))) stop("Parameter vector must be named!")
  
  # Convert parameter vector to a list for simnetfromtap:
  if (tmatch_type_pem != "no"){
    pems_low <- lapply(networks, function(x) x$pems[[1]])
    pems_high <- lapply(networks, function(x) x$pems[[2]])
    npems_low <- length(unique(unlist(lapply(pems_low, colnames))))
    npems_high <- length(unique(unlist(lapply(pems_high, colnames))))
    parList <- param_vec2list(params, n = npems_low, m = npems_high, fit.delta=fit.delta)
  } else {
    parList <- param_vec2list(params, n=0, m=0, fit.delta=fit.delta)
  }
  
  
  # Compute objective function value for each network:
  obj <- numeric(length(networks))
  for (i in 1:length(networks)) {
    paramsList <- parList
    # keep only latents that are used in the network:
    if (tmatch_type_pem != "no"){
      paramsList$lat_low <- parList$lat_low[which(names(parList$lat_low) %in% colnames(networks[[i]]$pems$low))]
      paramsList$lat_high <- parList$lat_high[which(names(parList$lat_high) %in% colnames(networks[[i]]$pems$high))]
    }
    I_mat <- simnetfromtap(traits = networks[[i]]$traits, abuns = networks[[i]]$abuns,
                           paramsList = paramsList, pems = networks[[i]]$pems,
                           tmatch_type_pem = tmatch_type_pem, tmatch_type_obs = tmatch_type_obs)
    if (!is.null(TmatchMatrixList)){ # check that there are no observed interactions where you have 0 probability!
      zeros1 <- which(TmatchMatrixList[[i]] == 0, arr.ind=FALSE)
      obs1 <- if (length(zeros1) > 0) networks[[i]]$web[zeros1] else 0
      if (any(obs1 > 0)) stop("The TmatchMatrixList contains zeros in places where observations were made. That leads to -Inf likelihoods.\n Please replace by a positive value substantially smaller than 1/prod(dim(web))!")
      I_mat <- I_mat * TmatchMatrixList[[i]]
      I_mat <- I_mat/sum(I_mat)
    }
    if (obj_function == "multinom") {
      obj[i] <- dmultinom(as.vector(networks[[i]]$web), size = sum(networks[[i]]$web),
                          prob = as.vector(I_mat), log = TRUE)
    } 
    
    if (obj_function == "least squares"){ # least squares:
      obj[i] <- sum((log(as.vector(networks[[i]]$web) + 1) - log(as.vector(I_mat) * sum(networks[[i]]$web) + 1))^2)
    }
    
    if (obj_function == "bjorn"){
      obj[i] <- suppressWarnings(-bjornloglik(web=networks[[i]]$web, P=I_mat))
    }
    
  } # end for loop
  
  if (obj_function == "multinom") {
    # add LASSO shrinkage (only for linear combination parameters of PEMs, *not* for trait matching parameters):
    sumOfAs <- if (tmatch_type_pem == "no") 0 else sum(abs(parList[[1]])) + sum(abs(parList[[2]]))
    out <- -sum(obj) + lambda * sumOfAs # note: - to use default minimisation to maximise likelihood
  } 
  
  if (obj_function == "least squares" | obj_function == "bjorn") out <- sum(obj) # note: no "-", as now the objective is the difference and hence to be minimised
  
  return(out)
}


#' @rdname internalFunctions
latent_cor <- function(true_pars, fitted_pars, pems_low, pems_high){
  # Compute correlation between true and fitted latent traits
  # for lower and higher trophic levels of a network
  
  # Select parameters
  true_lin_low <- true_pars$lat_low[which(names(true_pars$lat_low) %in% colnames(pems_low))]
  true_lin_high <- true_pars$lat_high[which(names(true_pars$lat_high) %in% colnames(pems_high))]
  fitted_lin_low <- fitted_pars$lat_low[which(names(fitted_pars$lat_low) %in% colnames(pems_low))]
  fitted_lin_high <- fitted_pars$lat_high[which(names(fitted_pars$lat_high) %in% colnames(pems_high))]
  true_shift <- true_pars$pem_shift
  fitted_shift <- fitted_pars$pem_shift
  
  # Calculate latent traits
  true_lat_low <- as.vector(scale(rowSums(matrix(true_lin_low, nrow = nrow(pems_low),
                                                 ncol = ncol(pems_low), byrow = TRUE) * pems_low)))
  true_lat_high <- as.vector(scale(rowSums(matrix(true_lin_high, nrow = nrow(pems_high),
                                                  ncol = ncol(pems_high), byrow = TRUE) * pems_high) + true_shift))
  fitted_lat_low <- as.vector(scale(rowSums(matrix(fitted_lin_low, nrow = nrow(pems_low),
                                                   ncol = ncol(pems_low), byrow = TRUE) * pems_low)))
  fitted_lat_high <- as.vector(scale(rowSums(matrix(fitted_lin_high, nrow = nrow(pems_high),
                                                    ncol = ncol(pems_high), byrow = TRUE) * pems_high) + fitted_shift))
  
  # Correlation between true and fitted traits
  out <- c("low" = cor(true_lat_low, fitted_lat_low), "high" = cor(true_lat_high, fitted_lat_high))
  return(out)
}


#' @rdname internalFunctions
web_indices <- function(web, web_dim, indices) {
  # Convert an interaction web vector back to a matrix and compute network indices
  web <- matrix(web, nrow = web_dim[1], ncol = web_dim[2])
  networklevel(web, index = indices)
}


#' @rdname internalFunctions
refit_params <- function(tapnet, fit, fitted_I_mat) {
  # Simulate new networks from a fitted tapnet object, re-fit on the simulated network
  # and output the parameter values
  for (i in 1:length(tapnet$networks)) {
    web <- tapnet$networks[[i]]$web
    tapnet$networks[[i]]$web <- matrix(rmultinom(1, sum(web), fitted_I_mat[[i]]), nrow = nrow(web),
                                       ncol = ncol(web))
  }
  refit <- fit_tapnet(tapnet, tmatch_type_pem = fit$tmatch_type_pem, tmatch_type_obs = fit$tmatch_type_obs,
                      lambda = fit$lambda, method = fit$method, maxit = fit$maxit)
  par_select <- refit$par_opt[3:5]
  names(par_select) <- NULL
  return(unlist(par_select))
}
