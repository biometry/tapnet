#' Predict from fitted tapnet to new data
#'
#' Function allows direct use of data prepared for tapnet analysis by other statistical methods, e.g. regression approaches
#'
#' The fitted tapnet object contains the estimated parameters, describing how traits, abundance and phylogeny play together to produce the network(s) used for fitting. This information is now used to predict interaction probabilities for a new network. Accordingly, we need to know this new network's species abundances (an input to the function), PEMs and traits. The latter are computed based on the information contained in the tapnet object (which is linked in by attribute reference). For new species, their PEMs are computed and the network is simulated, using \code{\link{simnetfromtap}}, for the information provided.
#'
#' @aliases predict_tapnet
#' 
#' @param fit results of applying fit_tapnet to the tapnet object;
#' @param abuns named list of two entries ("low" and "high"), containing a species-named vector of abundances of the new network
#' @param tapnet optional name of a tapnet object containing traits, phylogeny etc.; these are not stored in \option{fit}, but rather it is assumed that a tapnet object with the name stored in \option{fit} is available in the global environment. That may not be the case, e.g. when simulating networks. In this case, \option{tapnet} provides the required tapnet-object.
#' 
#' @return A matrix of predicted interaction probabilities, summing to 1. This would need to be multiplied by the total number of interactions in the new network to be comparable to the observations. 
#' 
#' @references Benadi et al. in prep
#'
#' @author Ruth Stephan, Gita Benadi and Carsten Dormann <carsten.dormann@biom.uni-freiburg.de>
#'
#' @examples
#' \dontrun{
#'   data(Tinoco)
#'   tap <- make_tapnet(tree_low = plant_tree, tree_high = humm_tree, networks = networks[1:2], 
#'          traits_low = plant_traits, traits_high = humm_traits, npems_lat = 4)
#'   fit <- fit_tapnet(tap) # uses two networks for fitting!
#'   gof_tapnet(fit)
#'   # predict to omitted forest network's abundances:
#'   pred1 <- predict_tapnet(fit, abuns=list("low"=plant_abun[[3]], "high"=humm_abun[[3]] )) 
#'   cor(as.vector(pred1*sum(networks[[3]])), as.vector(networks[[3]])) 
#' }
#' 
#' @export
predict_tapnet <- function(#tapnet, # A tapnet object upon which the prediction is based, created with simulate_tapnet
  fit, # Results of applying fit_tapnet to the tapnet object
  abuns,  # Abundances of lower and higher trophic level species in the predicted network (list of two vectors with named elements); if NULL, taken from tapnet object
  tapnet=NULL # optional: provide tapnet object if fit does not contain a valid link to an object (such as when the object was simulated but is not availabe in the global environment)
) {
  # Predict an interaction network with given species composition and abundances
  # based on parameter estimates from a fitted tapnet object
  
  if (is.null(tapnet)) tapnet <- get(attr(fit, "tapnet_name")) # uses attribute name to get the tapnet object
  
  # Check input format:
  if (class(tapnet) != "tapnet")
    stop("'tapnet' must be an object produced by function 'simulate_tapnet' or alike.")
  if (class(fit) != "fitted.tapnet")
    stop("'fit' must be an object produced by function 'fit_tapnet'.")
  # TODO: Check that "tapnet" is the object that was used to produce "fit".
  
  if (class(abuns) != "list" | length(abuns) != 2)
    stop("'abuns' must be a list of lower and higher trophic level abundances.")
  
  IsNamedVector <- function(vec) {
    is.vector(vec) & is.numeric(vec) & !is.null(names(vec)) & !any(is.na(names(vec)))
  } 
  
  if (!IsNamedVector(abuns[[1]]) | !IsNamedVector(abuns[[2]]))
    stop("Abundances must be vectors with named elements (species names).")
  
  # Check that all species occur in the phylogenetic tree:
  if (!all(names(abuns[[1]]) %in% tapnet$trees$low$tip.label) |
      !all(names(abuns[[2]]) %in% tapnet$trees$high$tip.label))
    stop("All species in the abundance list must appear in their respective phylogenetic tree.
			 Please make sure that species names match.")
  
  # Check that observed traits are available for all species:
  if (!is.null(tapnet$traits_all$low) & !is.null(tapnet$traits_all$high)) {
    if (!all(names(abuns[[1]]) %in% rownames(tapnet$traits_all$low)) |
        !all(names(abuns[[2]]) %in% rownames(tapnet$traits_all$high)))
      stop("Traits for some species are missing. Please make sure all species in the abundance
				 list have corresponding trait values. If necessary, add trait values to 'traits_all' in
				 the 'tapnet' object.")
  }
  
  # Calculate PEMs:
  pems <- list()
  pems_all_low <- pems_from_tree(tapnet$trees$low)
  pems$low  <- pems_all_low[rownames(pems_all_low) %in% names(abuns[[1]]), ]
  pems_all_high <- pems_from_tree(tapnet$trees$high)
  pems$high <- pems_all_high[rownames(pems_all_high) %in% names(abuns[[2]]), ]
  
  
  # Select those PEMs for which parameters were fitted:
  paramsList <- fit$par_opt
  pems$low <-  pems$low[which(names(pems$low) %in% names(paramsList$lat_low))]
  pems$high <-  pems$high[which(names(pems$high) %in% names(paramsList$lat_high))]
  
  # sorted alphabetically, because input is!
  
  
  # Select observed traits:
  traits <- list()
  traits_all_low <- tapnet$traits_all$low
  traits_all_high <- tapnet$traits_all$high
  if (!is.null(traits_all_low)){ # only if there are observed traits at all!
    traits$low <-  traits_all_low[rownames(traits_all_low) %in% names(abuns[[1]]), , drop = FALSE]
    traits$high <-  traits_all_high[rownames(traits_all_high) %in% names(abuns[[2]]), , drop = FALSE]
    # sort alphabetically!
    traits$low <- traits$low[order(rownames(traits$low)),, drop=FALSE]
    traits$high <- traits$high[order(rownames(traits$high)), , drop=FALSE]  
  }
  
  # Predict the matrix of interaction probabilities:
  predicted_Imat <- simnetfromtap(traits = traits, 
                                  abuns = abuns, 
                                  paramsList = paramsList, 
                                  pems = pems, 
                                  tmatch_type_pem = fit$tmatch_type_pem, 
                                  tmatch_type_obs = fit$tmatch_type_obs)
  
  return(predicted_Imat)
}
