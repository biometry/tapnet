#' Creates a simulated along with all parameters, abundances, traits and phylogeny used
#'
#' Simulation function to produce a tapnet object, i.e. one or more networks along with their descriptors, abundances, traits, phylogenetic information
#'
#' This function was written to explore the fitting of networks by different methods. Hence, we first have to simulate such networks. It is a very nice starting point for simulations, but irrelevant for tapnet-based analyses. The function internally sets up all the parameters and then calls \code{\link{simnetfromtap}} for the simulation of the actual network. Lots of options means, regrettably, lots of decisions for the user.
#'
#' @aliases simulate_tapnet
#' 
#' @param nlower species number in lower trophic level;
#' @param nhigher species number in higher trophic level;
#' @param ntraits_nopem number of phylogenetically uncorrelated traits (for each level);
#' @param ntraits_pem number of phylogenetically correlated traits (for each level);
#' @param pem_noise noise (sd of normal dist) to put on PEMs;
#' @param abuns abundances set to "lognormal", "equal" or a list of two abundance vectors;
#' @param meanlog parameters of the log-normal distribution for drawing abundances;
#' @param sdlog same as before, but width;
#' @param tmatch_type_pem  type of trait matching function for latent traits;
#' @param tmatch_type_obs  type of trait matching function for observed traits (can be a vector as long as there are traits);
#' @param npems_lat number of phylogenetic eigenvectors to be used to construct latent traits. If NULL, all eigenvectors will be used;
#' @param lat_low  vector of PEM linear combination parameters for lower trophic level; if 1, all values will be set to 1, if "random", values will be drawn from a uniform dist;
#' @param lat_high  same for higher trophic level;
#' @param tmatch_width_pem width of trait matching function for latent (PEM-based) traits;
#' @param pem_shift  shift parameter for latent trait matching;
#' @param tmatch_width_obs width of trait matching function for observed traits, can be a single value or a vector of length ntraits_nopem + ntraits_pem;
#' @param Nwebs number of webs to be simulated;
#' @param prop_species  proportion of species in the phylogeny that appear in each web (species are drawn randomly). With multiple networks, this allows to create networks with partly overlapping species composition.
#' @param new_abuns  If abuns = "lognormal" and Nwebs > 1, should new abundances be drawn for each web?
#' @param Nobs number of observed interactions per web
#' 
#' @return A tapnet object, with a highly nested structure. There are six entries at the top level: trees, traits_all, networks, sim_params and the two tmatch types. Within networks, there is a list of 5 entries (for each of the \option{Nweb} networks): abundances, traits, PEMs, web and I_mat. "I_mat" is the actual output from simnetfromtap, while "web" is a single draw from a multinomial distribution with I_mat as probabilities and \option{Nobs} as size.
#' 
#' 
#' 
#' @references Benadi et al. in prep
#'
#' @author Gita Benadi <gita.benadi@biom.uni-freiburg.de>, Jochen Fründ <jochen.fruend@biom.uni-freiburg.de> and Carsten Dormann <carsten.dormann@biom.uni-freiburg.de>
#'
#' @examples
#' tapnet <- simulate_tapnet(nlower=10, nhigher=20, ntraits_nopem=2, ntraits_pem=0) 
#' # a minimal call of simulate_tapnet
#' str(tapnet, 1) # the structure at the first level
#' str(tapnet, 2) # the structure at the first and second level
#'
#' @export
simulate_tapnet <- function(nlower, # Species number in lower trophic level
                            nhigher, # Species number in higher trophic level
                            ntraits_nopem, # Number of phylogenetically uncorrelated traits (for each level)
                            ntraits_pem,  # Number of phylogenetically correlated traits (for each level)
                            pem_noise = 0.5, # Noise (sd of normal dist) to put on PEMs
                            abuns = "lognormal", # set abundances to "lognormal", "equal" or a list of two abundance vectors
                            meanlog = 0, # Parameters of the log-normal distribution for drawing abundances
                            sdlog = 1,
                            tmatch_type_pem = "normal", # type of trait matching function for latent traits
                            tmatch_type_obs = "normal",  # type of trait matching function for observed traits
                            # (can be a vector)
                            npems_lat = NULL, # Number of phylogenetic eigenvectors to be used to construct
                            # latent traits. If NULL, all eigenvectors will be used.
                            lat_low = 1, # Vector of PEM linear combination parameters for lower trophic level
                            # if 1, all values will be set to 1, if "random", values will be drawn from a uniform dist.
                            lat_high = 1, # Same for higher trophic level
                            tmatch_width_pem = 1, # Width of trait matching function for latent (PEM-based) traits
                            pem_shift = 0, # Shift parameter for latent trait matching
                            tmatch_width_obs = 1, # Width of trait matching function for observed traits
                            # can be a single value or a vector of length ntraits_nopem + ntraits_pem
                            Nwebs = 1, # Number of webs to be simulated
                            prop_species = 1, #  Proportion of species in the phylogeny
                            # that appear in each web (species are drawn randomly). With multiple networks,
                            # this allows to create networks with partly overlapping species composition.
                            new_abuns = FALSE, # If abuns = "lognormal" and Nwebs > 1, should new abundances be drawn
                            # for each web?
                            Nobs = 1000 # Number of observed interactions per web
) {
  # Simulate traits, abundances, phylogenies and interaction networks
  
  # TODO: Remove ntraits_nopem? Phylogenetically uncorrelated traits can also be created
  # by putting a lot of noise on the PEMs.
  
  # TODO: Scale and center simulated trait values, to make trait matching parameters more comparable?
  
  # 1. Phylogenies
  
  # Simulate phylogenies:
  tree_low <- pbtree(n = nlower, scale = 1, nsim = 1)
  tree_low$tip.label <- paste0("L", 1:nlower)
  tree_high <- pbtree(n = nhigher, scale = 1, nsim = 1)
  tree_high$tip.label <- paste0("H", 1:nhigher)
  
  # Calculate PEMs for all species:
  pems_low <- pems_from_tree(tree_low)
  pems_high <- pems_from_tree(tree_high)
  
  # 2. Traits
  
  if (!(ntraits_nopem > 0 | ntraits_pem > 0)) {
    # Set trait matrices to NULL if there are no observed traits
    traits_low <- NULL
    traits_high <- NULL
  } else {
    # Check that number of PEM-based traits is not higher than number of PEMs:
    if (ntraits_pem > min(nlower, nhigher) - 1)
      stop("'ntraits_pem' must be smaller than min(nlower, nhigher)!")
    
    # Draw phylogenetically uncorrelated traits (from lognormal distribution):
    if (ntraits_nopem > 0) {
      traits_low <- matrix(rlnorm(nlower * ntraits_nopem), nrow = nlower, ncol = ntraits_nopem)
      traits_high <- matrix(rlnorm(nhigher * ntraits_nopem), nrow = nhigher, ncol = ntraits_nopem)
    }
    
    # Put noise on selected PEMs to create phylogenetically correlated traits:
    if (ntraits_pem > 0) {
      traits_low_pem <- pems_low[, 1:ntraits_pem] + rnorm(nlower * ntraits_pem,
                                                          mean = 0, sd = pem_noise)
      traits_high_pem <- pems_high[, 1:ntraits_pem] + rnorm(nhigher * ntraits_pem,
                                                            mean = 0, sd = pem_noise)
      if (ntraits_nopem > 0) {
        traits_low <- cbind(traits_low, traits_low_pem)
        traits_high <- cbind(traits_high, traits_high_pem)
      } else {
        traits_low <- as.matrix(traits_low_pem)
        traits_high <- as.matrix(traits_high_pem)
      }
    }
    # Add species names
    rownames(traits_low) <- tree_low$tip.label
    rownames(traits_high) <- tree_high$tip.label
    # Number traits
    colnames(traits_low) <- paste0("T", 1:ncol(traits_low))
    colnames(traits_high) <- paste0("T", 1:ncol(traits_high))
  }
  
  
  # 3. Species composition, PEM and trait selection
  # (for webs with less than the full set of species, i.e. prop_species < 1)
  
  if (prop_species <= 0 | prop_species > 1)
    stop("'prop_species' must be a proportion (i.e. between 0 and 1, excluding 0).")
  
  if (prop_species < 1) {
    nspec_low <- round(prop_species * nlower)
    nspec_high <- round(prop_species * nhigher)
    if (nspec_low < 2 | nspec_high < 2) stop("The chosen proportion of species per network is too low.")
    spec_low <- replicate(Nwebs, sample(tree_low$tip.label, nspec_low, replace = F), simplify = F)
    spec_low <- lapply(spec_low, function(x) x[order(nchar(x), x)])
    spec_high <- replicate(Nwebs, sample(tree_high$tip.label, nspec_high, replace = F), simplify = F)
    spec_high <- lapply(spec_high, function(x) x[order(nchar(x), x)])
    pems_low <- lapply(spec_low, select_relevant_pems, tree = tree_low)
    pems_high <- lapply(spec_high, select_relevant_pems, tree = tree_high)
    if (!is.null(traits_low)) {
      trait_list_low <- list()
      trait_list_high <- list()
      for (i in 1:Nwebs) {
        trait_list_low[[i]] <- traits_low[which(rownames(traits_low) %in% spec_low[[i]]), , drop = F]
        trait_list_high[[i]] <- traits_high[which(rownames(traits_high) %in% spec_high[[i]]), , drop = F]
      }
    }
  } else {
    nspec_low <- nlower
    nspec_high <- nhigher
    if (!is.null(traits_low)) {
      trait_list_low <- rep(list(traits_low), Nwebs)
      trait_list_high <- rep(list(traits_high), Nwebs)
    }
    pems_low <- rep(list(pems_low), Nwebs)
    pems_high <- rep(list(pems_high), Nwebs)
  }
  if (!is.null(npems_lat)) {
    pems_low <- lapply(pems_low, function(x) x[, 1:npems_lat, drop = F])
    pems_high <- lapply(pems_high, function(x) x[, 1:npems_lat, drop = F])
  }
  pems_all_low <- unique(unlist(lapply(pems_low, names)))
  pems_all_low <- pems_all_low[order(nchar(pems_all_low), pems_all_low)]
  pems_all_high <- unique(unlist(lapply(pems_high, names)))
  pems_all_high <- pems_all_high[order(nchar(pems_all_high), pems_all_high)]
  
  # 4. Abundances
  
  if (!(abuns %in% c("equal", "lognormal", "log-normal") | is.list(abuns)))
    stop("'abuns' must be a list of abundances (low, high), \"equal\", or \"lognormal\".")
  
  if (abuns == "equal" | is.list(abuns) | (abuns %in% c("lognormal", "log-normal") & new_abuns == F)) {
    if (abuns == "equal") {
      abun_low <- rep(1 / nspec_low, nspec_low)
      abun_high <- rep(1 / nspec_high, nspec_high)
    }
    if (is.list(abuns)) {
      if (length(abuns[[1]]) != nspec_low | length(abuns[[2]]) != nspec_high)
        stop("Abundance vectors don't match number of species.")
      abun_low <- abuns[[1]]
      abun_high <- abuns[[2]]
    }
    if (abuns %in% c("lognormal", "log-normal")) {
      abun_low <- rlnorm(nspec_low, meanlog = meanlog, sdlog = sdlog)
      abun_low <- abun_low / sum(abun_low)
      abun_high <- rlnorm(nspec_high, meanlog = meanlog, sdlog = sdlog)
      abun_high <- abun_high / sum(abun_high)
    }
    abun_low <- rep(list(abun_low), Nwebs)
    abun_high <- rep(list(abun_high), Nwebs)
  } else {
    abun_low <- replicate(Nwebs, rlnorm(nspec_low, meanlog = meanlog, sdlog = sdlog), simplify = F)
    abun_low <- lapply(abun_low, function(x) x / sum(x))
    abun_high <- replicate(Nwebs, rlnorm(nspec_high, meanlog = meanlog, sdlog = sdlog), simplify = F)
    abun_high <- lapply(abun_high, function(x) x / sum(x))
  }
  # Add species names
  if (prop_species < 1) {
    for (i in 1:Nwebs) {
      names(abun_low[[i]]) <- spec_low[[i]]
      names(abun_high[[i]]) <- spec_high[[i]]
    }
  } else {
    for (i in 1:Nwebs) {
      names(abun_low[[i]]) <- tree_low$tip.label
      names(abun_high[[i]]) <- tree_high$tip.label
    }
  }
  
  
  # 5. Parameters
  
  # Latent trait combination parameters
  nlat_low <- length(pems_all_low)
  nlat_high <- length(pems_all_high)
  
  if (length(lat_low) == 1) {
    if (!(lat_low == 1 | lat_low == "random"))
      stop("'lat_low' must be either 1, 'random' or a numeric vector.")
  } else {
    if (!is.numeric(lat_low))
      stop("'lat_low' must be either 1, 'random' or a numeric vector.")
  }
  if (length(lat_high) == 1) {
    if (!(lat_high == 1 | lat_high == "random"))
      stop("'lat_high' must be either 1, 'random' or a numeric vector.")
  } else {
    if (!is.numeric(lat_high))
      stop("'lat_high' must be either 1, 'random' or a numeric vector.")
  }
  
  if (length(lat_low) == 1) {
    if (lat_low == 1) {
      lat_low <- rep(1, nlat_low)
    } else {
      if (lat_low == "random") lat_low <- runif(nlat_low)
    }
  } else {
    if (prop_species < 1)
      stop("Values of 'lat_low' and 'lat_high' can only be specified when 'prop_species = 1'")
    if (length(lat_low) != nlat_low) stop("'lat_low' must be of length nlower - 1.")
  }
  
  if (length(lat_high) == 1) {
    if (lat_high == 1) {
      lat_high <- rep(1, nlat_high)
    } else {
      if (lat_high == "random") lat_high <- runif(nlat_high)
    }
  } else {
    if (prop_species < 1)
      stop("Values of 'lat_low' and 'lat_high' can only be specified when 'prop_species = 1'")
    if (length(lat_high) != nlat_high) stop("'lat_high' must be of length nhigher - 1.")
  }
  
  names(lat_low) <- pems_all_low
  names(lat_high) <- pems_all_high
  
  # Trait matching parameters
  if (tmatch_width_pem <= 0 | any(tmatch_width_obs <= 0)) stop("'tmatch_width' parameters must be positive.")
  if (length(tmatch_width_obs) > 1 & length(tmatch_width_obs) != ntraits_nopem + ntraits_pem) {
    stop("'tmatch_width_obs' must either be a single value or a vector of length ntraits_nopem + ntraits_pem.")
  }
  if ((length(tmatch_width_obs) == 1) & (ntraits_nopem + ntraits_pem > 1)) {
    tmatch_width_obs <- rep(tmatch_width_obs, ntraits_nopem + ntraits_pem)
  }
  
  paramsList <- list(lat_low = lat_low, lat_high = lat_high, pem_shift = pem_shift,
                     tmatch_width_pem = tmatch_width_pem, tmatch_width_obs = tmatch_width_obs)
  
  # 6. Simulate network(s)
  
  # If abundances and/or species composition vary between networks,
  # draw separate interaction probability matrices, otherwise draw just one
  webs <- list()
  for (i in 1:Nwebs) {
    webs[[i]] <- list()
    webs[[i]]$abuns <- list(low = abun_low[[i]], high = abun_high[[i]])
    if (is.null(traits_low)) {
      webs[[i]]$traits <- list(low = traits_low, high = traits_high)
    } else {
      webs[[i]]$traits <- list(low = trait_list_low[[i]], high = trait_list_high[[i]])
    }
    webs[[i]]$pems <- list(low = pems_low[[i]], high = pems_high[[i]])
    if (prop_species < 1) {
      paramsList$lat_low <- lat_low[which(names(lat_low) %in% colnames(pems_low[[i]]))]
      paramsList$lat_high <- lat_high[which(names(lat_high) %in% colnames(pems_high[[i]]))]
    }
    webs[[i]]$I_mat <- simnetfromtap(webs[[i]]$traits, webs[[i]]$abuns, paramsList, webs[[i]]$pems,
                                     tmatch_type_pem, tmatch_type_obs)
    webs[[i]]$web <- matrix(rmultinom(1, Nobs, webs[[i]]$I_mat), nrow = nrow(webs[[i]]$I_mat),
                            ncol = ncol(webs[[i]]$I_mat))
    dimnames(webs[[i]]$web) <- dimnames(webs[[i]]$I_mat)
  }
  
  
  # 7. Combine everything into a tapnet object
  
  paramsList$lat_low <- lat_low
  paramsList$lat_high <- lat_high
  
  sim_tapnet <- list(trees = list(low = tree_low, high = tree_high),
                     traits_all = list(low = traits_low, high = traits_high), networks = webs,
                     sim_params = paramsList, tmatch_type_pem = tmatch_type_pem,
                     tmatch_type_obs = tmatch_type_obs)
  class(sim_tapnet) <- "tapnet"
  return(sim_tapnet)
}