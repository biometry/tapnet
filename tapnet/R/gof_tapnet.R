#' Goodness-of-fit of a tapnet fit
#'
#' Provides various measures to describe how well the tapnet model fits the data
#'
#' This is a function particularly interesting for simulated data, when the true parameters are known. In this case, GOF for fitted latent traits and so forth are computed.
#'
#' @aliases gof_tapnet
#' 
#' @param fit results of applying fit_tapnet to the tapnet object;
#' @param tapnet a tapnet object created with simulate_tapnet or make_tapnet; if NULL, the name stored in the attributes of fit is used to access an object of that name in the global environment; can be used e.g. in simulations, when the tapnet object is renamed relative to the one fitted;
#' @param indices network indices to compare between observed and fitted network; calls \code{\link[bipartite]{networklevel}};
#' @param nrep Number of networks to simulate for indices comparison; these are draws from the fitted multinomial distribution;
#' @param se_refit logical; should standard errors for the parameters be calculated using parametric bootstrap (refitting on simulated data)? Defaults to FALSE because it's very slow (i.e. takes hours).
#' @param se_nsim number of simulations for parametric bootstrap (ignored unless the previous argument is set to TRUE).
#' 
#' @return A list of goodness-of-fit measures: bc_sim_web are the Bray-Curtis similarities between fitted and observed network; cor_web are Spearman correlations between fitted and observed; and net_indices compute the selected network indices for fitted and observed networks. If more than one network is used for fitting, all these measures are returned for all networks (as vector or list under the respective label). See example.
#' 
#' 
#' @references Benadi et al. in prep
#'
#' @author Gita Benadi <gita.benadi@biom.uni-freiburg.de> and Carsten Dormann <carsten.dormann@biom.uni-freiburg.de>
#'
#' @examples
#' \dontrun{
#'   data(Tinoco)
#'   tap <- make_tapnet(tree_low = plant_tree, tree_high = humm_tree, networks = networks[2:3], traits_low = plant_traits, traits_high = humm_traits, npems_lat = 4)
#'   fit <- fit_tapnet(tap) # uses networks 2 and 3 for fitting!
#'   gof_tapnet(fit)
#' }
#'
#' @export
gof_tapnet <- function(
  fit, # Results of applying fit_tapnet to the tapnet object
  tapnet=NULL, # A tapnet object created with simulate_tapnet or make_tapnet; if NULL, the name stored in the attributes of fit is used to access an object of that name in the global environment; can be used e.g. in simulations, when the tapnet object is renamed relative to the one fitted
  indices = c("connectance", "NODF", "weighted NODF", "H2"), # Network indices to compare
  nrep = 1000, # Number of networks to simulate for indices comparison
  se_refit = FALSE, # Should standard errors for the parameters be calculated using
  # parametric bootstrap (refitting on simulated data)?
  # Defaults to FALSE because it's very slow (i.e. takes *hours*).
  se_nsim = 1000 # Number of simulations for parametric bootstrap
) {
  # Assess goodness-of-fit in different ways
  if (is.null(tapnet)) tapnet <- get(attr(fit, "tapnet_name")) # uses attribute name to get the tapnet object
  
  gof <- list()
  
  # Check input format
  if (class(tapnet) != "tapnet")
    stop("'tapnet' must be an object produced by function 'simulate_tapnet'.")
  if (class(fit) != "fitted.tapnet")
    stop("'fit' must be an object produced by function 'fit_tapnet'.")
  # TODO: Check that "tapnet" is the object that was used to produce "fit".
  
  
  # 1. Compare fitted and true interaction matrices
  
  nwebs <- length(tapnet$networks)
  fitted_I_mat <- list()
  fitted_webs <- list()
  gof$bc_sim_web <- numeric(nwebs)
  gof$cor_web <- numeric(nwebs)
  gof$net_indices <- list()
  
  for (i in 1:nwebs) {
    web <- tapnet$networks[[i]]
    paramsList <- fit$par_opt
    paramsList$lat_low <- paramsList$lat_low[which(names(paramsList$lat_low) %in% colnames(web$pems$low))]
    paramsList$lat_high <- paramsList$lat_high[which(names(paramsList$lat_high) %in% colnames(web$pems$high))]
    fitted_I_mat[[i]] <- simnetfromtap(traits = web$traits, abuns = web$abuns, params = paramsList,
                                       pems = web$pems, tmatch_type_pem = fit$tmatch_type_pem,
                                       tmatch_type_obs = fit$tmatch_type_obs)
    # a) Compare entries of fitted and true interaction matrices
    fitted_webs[[i]] <- fitted_I_mat[[i]] * sum(web$web)
    gof$bc_sim_web[i] <- 1 - vegdist(matrix(c(as.vector(web$web), as.vector(fitted_webs[[i]])),
                                            nrow = 2, byrow = TRUE), method = "bray")[1]
    gof$cor_web[i] <- cor(as.vector(web$web), as.vector(fitted_webs[[i]]), method = "spearman") 
    # b) Compare network metrics of fitted and true webs
    true_web <- tapnet$networks[[i]]$web
    fitted_sim <- rmultinom(nrep, sum(true_web), fitted_I_mat[[i]]) # repeatedly sample from fitted multinomial distribution
    fitted_indices <- t(apply(fitted_sim, 2, web_indices, web_dim = dim(true_web), indices = indices))
    obs <- networklevel(true_web, index = indices)
    gof$net_indices[[i]] <- data.frame(Index = colnames(fitted_indices),
                                       Observed = obs,
                                       Mean = apply(fitted_indices, 2, mean, na.rm = TRUE),
                                       Median = apply(fitted_indices, 2, median, na.rm = TRUE),
                                       q2.5 = apply(fitted_indices, 2, quantile, probs = 0.025, na.rm = TRUE),
                                       q97.5 = apply(fitted_indices, 2, quantile, probs = 0.975, na.rm = TRUE))
    rownames(gof$net_indices[[i]]) <- NULL
  }
  
  # 2. Compare fitted and true parameters, estimate uncertainty of fitted parameters
  
  # Correlation of true and fitted latent traits (only for simulated networks)
  if ("sim_params" %in% names(tapnet)) { # Check that this is a simulated network
    pem_list_low <- lapply(tapnet$networks, function(x) x$pems[[1]])
    pem_list_high <- lapply(tapnet$networks, function(x) x$pems[[2]])
    gof$latent_cor <- mapply(latent_cor, pem_list_low, pem_list_high,
                             MoreArgs = list(true_pars = tapnet$sim_params, fitted_pars = fit$par_opt))
  }
  
  # Parameter standard errors calculated using the Hessian matrix
  # (Not for latent trait coefficients)
  if (!is.null(fit$opt$hessian)) {
    gof$se_hess <- suppressWarnings(sqrt(diag(solve(-fit$opt$hessian)))[-(1:(length(paramsList$lat_low) +
                                                                               length(paramsList$lat_high)))])
  }
  
  # Parameter standard errors calculated using parametric bootstrap
  if (se_refit == TRUE) {
    refit <- t(replicate(se_nsim, refit_params(tapnet, fit, fitted_I_mat)))
    gof$se_refit <- apply(refit, 2, sd)
  }
  
  return(gof)
  
}
