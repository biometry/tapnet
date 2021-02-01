#' Helper functions for tapnet
#'
#' Lower-level, non-exported functions to be called by the main tapnet functions
#'
#' They do roughly the following: 
#' \describe{
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
#' @param lambda LASSO shrinkage parameter to avoid collinearity issues when observed traits are phylogenetically correlated;
#' @param obj_function objective function, either "multinom" or anything else (currently anything else leads to OLS fitting);
#' @param true_pars parameters used for simulating the network;
#' @param fitted_pars parameters estimated by \code{\link{fit_tapnet}};
#' @param pems_low phylogenetic eigenvectors for the lower trophic level;
#' @param pems_high phylogenetic eigenvectors for the higher trophic level;
#' @param web data for an interaction network in vector form, e.g. from predict_tapnet;
#' @param web_dim vector of two numbers indicating the rows and column of the network matrix;
#' @param indices vector of names of network indices to compute; see \code{\link[bipartite]{networklevel}} for what is available;
#' @param tapnet a tapnet object;
#' @param fit a fitted tapnet;
#' @param fitted_I_mat the fitted I-matrix of a fitted tapnet object (I think).

#'
#' @references Benadi et al. in prep
#'
#' @author Gita Benadi <gita.benadi@biom.uni-freiburg.de>, Carsten Dormann <carsten.dormann@biom.uni-freiburg.de> and Jochen Fründ <jochen.fruend@biom.uni-freiburg.de>
#'
#'
internalFunctions <- function(){}# just here to display the name for the help page through Roxygen correctly; grrrr

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
                   err = 1E-5 # "baseline" probability of match, even if traits do not match at all
){# Calculate interaction probabilities based on trait matching
  # shift
  delta_t <- delta_t + shift
  
  # lognormal distribution with mode shifted to zero:
  if(type == "shiftlnorm") out <- dlnorm(delta_t + width, meanlog = log(width) + 1) + err
  
  # normal distribution:
  if(type == "normal") out <- dnorm(delta_t, mean = 0, sd = width) + err
  
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
  par_list$lat_low <- params[1:n] # PEM linear combination parameters (lower level)
  par_list$lat_low[1] <- exp(params[1]) # Force positive value
  par_list$lat_high <- params[(n + 1):(n + m)] # PEM linear combination parameters (higher level)
  par_list$pem_shift <- params[n + m + 1] # PEM shift
  par_list$tmatch_width_pem <- exp(params[n + m + 2]) # PEM trait matching
  if (fit.delta){ # if delta is to be fitted, it is the last of params
    if (length(params) > (n + m + 2)) {
      par_list$tmatch_width_obs <- exp(params[(n + m + 3):(length(params)-1)]) # observed traits matching
    }
    par_list$delta <- plogis(params[length(params)])
  } else { # if delta is NOT fitted, all last parameters go to traits
    if (length(params) > (n + m + 2)) {
      par_list$tmatch_width_obs <- exp(params[(n + m + 3):(length(params))]) # observed traits matching
    }
  }
  return(par_list)
}

#' @rdname internalFunctions
loglik_tapnet <- function(params, # Parameters (a *named* vector)
                          networks, # the "networks" part of a tapnet object
                          tmatch_type_pem, # Type of trait matching function for latent traits
                          tmatch_type_obs, # Type(s) of trait matching functions for observed traits
                          lambda=0, # LASSO shrinkage parameter
                          obj_function = "multinom", # Objective function:
                          fit.delta=TRUE # either "multinom" (negative log-likelihood based on a multinomial distribution) or
                          # "sq_diff" (sum of squared differences between observed and predicted log(interactions + 1))
) {
  # Compute value of the objective function for fitting given the parameter values and "tap" data
  
  if (is.null(names(params))) stop("Parameter vector must be named!")
  
  # Convert parameter vector to a list for simnetfromtap:
  pems_low <- lapply(networks, function(x) x$pems[[1]])
  pems_high <- lapply(networks, function(x) x$pems[[2]])
  npems_low <- length(unique(unlist(lapply(pems_low, colnames))))
  npems_high <- length(unique(unlist(lapply(pems_high, colnames))))
  parList <- param_vec2list(params, n = npems_low, m = npems_high, fit.delta=fit.delta)
  
  # Compute objective function value for each network:
  obj <- numeric(length(networks))
  for (i in 1:length(networks)) {
    paramsList <- parList
    # keep only latents that are used in the network:
    paramsList$lat_low <- parList$lat_low[which(names(parList$lat_low) %in% colnames(networks[[i]]$pems$low))]
    paramsList$lat_high <- parList$lat_high[which(names(parList$lat_high) %in% colnames(networks[[i]]$pems$high))]
    I_mat <- simnetfromtap(traits = networks[[i]]$traits, abuns = networks[[i]]$abuns,
                           paramsList = paramsList, pems = networks[[i]]$pems,
                           tmatch_type_pem = tmatch_type_pem, tmatch_type_obs = tmatch_type_obs)
    if (obj_function == "multinom") {
      obj[i] <- dmultinom(as.vector(networks[[i]]$web), size = sum(networks[[i]]$web),
                          prob = as.vector(I_mat), log = TRUE)
    } else { # least squares:
      obj[i] <- sum((log(as.vector(networks[[i]]$web) + 1) - log(as.vector(I_mat) * sum(networks[[i]]$web) + 1))^2)
    }
  }
  if (obj_function == "multinom") {
    # add LASSO shrinkage (only for linear combination parameters of PEMs, *not* for trait matching parameters):
    sumOfAs <- sum(abs(parList[[1]])) + sum(abs(parList[[2]]))
    return(-sum(obj) + lambda * sumOfAs) # note: - to use default minimisation to maximise likelihood
  } else {
    return(sum(obj)) # note: no "-", as now the objective is the difference and hence to be minimised
  }
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
