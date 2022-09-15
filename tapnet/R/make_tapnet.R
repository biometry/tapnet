#' Constructs an object of type "tapnet"
#'
#' Collates networks, traits and phylogenies into a consistent data structure used for all other tapnet functions
#'
#' Tapnet objects are the starting point for almost all other tapnet functions. They contain the information on the species and the (quantitative) interaction network data.
#'
#' @aliases make_tapnet
#' @param tree_low phylogenetic tree of lower trophic level; required;
#' @param tree_high	 phylogenetic tree of higher trophic level; required;
#' @param networks a single or list of interaction network (as matrix); required; 
#' @param abun_low named abundance vector(s) for lower trophic level (single vector or list of vectors); optional;
#' @param abun_high named abundance vector(s) for higher trophic level; optional;
#' @param traits_high higher trophic level traits (species x traits matrix with row and column names; optional;
#' @param traits_low lower trophic level traits (species x traits matrix with row and column names); optional;
#' @param npems_lat number of phylogenetic eigenvectors to be used to construct latent traits. If NULL, all eigenvectors will be used.
#' @param use.all.pems option to force the function to use all phylogenetic eigenvectors, not only those useful for describing the specific network's species.
#'
#' @return A tapnet object, i.e. an thoroughly organised list with the inputs as entries. If multiple networks are provided, each has its own list entry, with PEMs, traits and abundances given for each network separately, in addition to the overall phylogenetic eigenvectors across all networks. See example for, well, for an example.
#'
#' @references Benadi et al. in prep
#'
#' @author Gita Benadi <gita.benadi@biom.uni-freiburg.de>
#'
#' @examples
#' data(Tinoco)
#' tapnet_web1 <- make_tapnet(tree_low = plant_tree, tree_high = humm_tree, networks = networks[1], 
#'                traits_low = plant_traits, traits_high = humm_traits, npems_lat = 4)
#' str(tapnet_web1) # show tapnet structure
#'
#' @export
make_tapnet <- function(tree_low, # phylogenetic tree of lower trophic level (required)
                        tree_high, # phylogenetic tree of higher trophic level (required)
                        networks, # a single interaction network (as matrix) or a list of networks (required)
                        abun_low = NULL, # named abundance vector(s) for lower trophic level
                        # (single vector or list of vectors, optional)
                        abun_high = NULL, # named abundance vector(s) for higher trophic level (optional)
                        traits_low = NULL, # lower trophic level traits (species x traits matrix with row and column names, optional)
                        traits_high = NULL, # higher trophic level traits (species x traits matrix with row and column names, optional)
                        npems_lat = NULL, # Number of phylogenetic eigenvectors to be used to construct latent traits. If NULL, all eigenvectors will be used.
                        use.all.pems=FALSE # sort through phylogenetic eigenvectors and use only those relevant for this network (see helper function select_relevant_pems).
) {
  # Construct an object of class "tapnet" from the supplied data
  
  # TODO: Maybe modify the function so it always performs all checks, collect messages
  # about problems and only then stops and outputs all error messages.
  
  # Check the input:
  
  # Phylogenies
  if (! (("phylo" %in% class(tree_low)) | "phylo" %in% class(tree_high) )) {
    stop("'tree_low' and 'tree_high' must be of class 'phylo'.
			 See the documentation of package 'ape' on how to convert phylogenetic trees to this format.")
  }
  
  # Networks
  if (!is.list(networks)) networks <- list(networks)
  mat_check <- sapply(networks, is.matrix)
  if (any(mat_check == FALSE)) stop("All elements of 'networks' must be matrices.")
  
  # Are all species named?
  for (i in 1:length(networks)) {
    if (is.null(rownames(networks[[i]])) | is.null(colnames(networks[[i]]))) {
      stop("All interaction networks must have species identities as row and column names.")
    }
  }
  
  spec_lower <- sort(unique(unlist(lapply(networks, rownames)))) # names of lower trophic level species in all networks
  spec_higher <- sort(unique(unlist(lapply(networks, colnames)))) # same for higher trophic level
  
  # Do all interacting species appear in the phylogenies?
  if (!all(spec_lower %in% tree_low$tip.label) | !all(spec_higher %in% tree_high$tip.label)) {
    stop("All species in the interaction networks must appear in their respective phylogeny.
			 Please make sure that species names match.")
  }
  
  # Sort species alphabetically
  for (i in 1:length(networks)) {
    networks[[i]] <- networks[[i]][order(rownames(networks[[i]])), order(colnames(networks[[i]]))]
  }
  
  # Abundances
  if (!is.null(abun_low)) {
    if (!is.list(abun_low)) abun_low <- list(abun_low)
    num_check <- sapply(abun_low, is.numeric)
    if (any(num_check == F)) stop("All elements of 'abun_low' must be numeric vectors.")
    # Do we have abundances for all networks?
    if (length(networks) != length(abun_low)) stop("Some networks have missing lower level abundances.")
    # Do we have abundances for all species? Sort species alphabetically.
    for (i in 1:length(networks)) {
      abun_low[[i]] <- abun_low[[i]][order(names(abun_low[[i]]))]
      if (!all(rownames(networks[[i]]) == names(abun_low[[i]]))) stop("Some lower trophic level species have missing abundances.")
    }
  } else {
    warning("No abundances for lower trophic level were provided. Using marginal totals instead.")
  }
  
  if (!is.null(abun_high)) {
    if (!is.list(abun_high)) abun_high <- list(abun_high)
    num_check <- sapply(abun_high, is.numeric)
    if (any(num_check == F)) stop("All elements of 'abun_high' must be numeric vectors.")
    # Do we have abundances for all networks?
    if (length(networks) != length(abun_high)) stop("Some networks have missing higher level abundances.")
    # Do we have abundances for all species? Sort species alphabetically.
    for (i in 1:length(networks)) {
      abun_high[[i]] <- abun_high[[i]][order(names(abun_high[[i]]))]
      if (!all(colnames(networks[[i]]) == names(abun_high[[i]]))) stop("Some higher trophic level species have missing abundances.")
    }
  } else {
    warning("No abundances for higher trophic level were provided. Using marginal totals instead.")
  }
  
  # TODO: Check that species names in abundance vectors and networks match and sort species alphabetically
  
  # Traits
  if (!is.null(traits_low) | !is.null(traits_high)) {
    if (is.null(traits_low) | is.null(traits_high)) {
      stop("Trait data must be supplied for both trophic levels.")
    }
    if (!is.matrix(traits_low) | !is.matrix(traits_high)) {
      stop("'traits_low' and 'traits_high' must be matrices.")
    }
    if (ncol(traits_low) != ncol(traits_high)) {
      stop("Numbers of traits of lower and higher trophic level must be equal.")
    }
    # Do all interacting species have trait data?
    if (!all(spec_lower %in% rownames(traits_low))) stop("Some lower level species have missing trait data.")
    if (!all(spec_higher %in% rownames(traits_high))) stop("Some higher level species have missing trait data.")
  }
  
  # Build tapnet object
  trees <- list(low = tree_low, high = tree_high)
  traits_all <- list(low = traits_low, high = traits_high)
  
  if (is.null(abun_low)) {
    abun_low <- lapply(networks, rowSums)
  }
  if (is.null(abun_high)) {
    abun_high <- lapply(networks, colSums)
  }
  
  webs <- list()
  for (i in 1:length(networks)) {
    webs[[i]] <- list()
    webs[[i]]$web <- networks[[i]]
    if (use.all.pems == FALSE){
      pems_web_low <- select_relevant_pems(tree_low, rownames(networks[[i]]))
    } else {
      pems_web_low <- pems_from_tree(tree_low)
    }
    pems_web_low <- pems_web_low[order(rownames(pems_web_low)),] # sort species alphabetically
    if (use.all.pems == FALSE){
      pems_web_high <- select_relevant_pems(tree_high, colnames(networks[[i]]))
    } else {
      pems_web_high <- pems_from_tree(tree_high)
    }
    pems_web_high <- pems_web_high[order(rownames(pems_web_high)),] # sort species alphabetically
    if (npems_lat == 0){ # use any PEM?
      warning("No phylogenetic information will be used (otherwise change option 'npems_lat').")
      pems_web_low <- NA
      pems_web_high <- NA
    } else {
        if (!is.null(npems_lat)) { # use stated PEMs
        if (npems_lat > ncol(pems_web_low)) {
          warning(paste("For web no.", i, ", only",  ncol(pems_web_low) , "PEMs will be used for the lower trophic level latent trait,\n since this network only has", nrow(webs[[i]]$web), "lower trophic level species."))
        } else {
          pems_web_low <- pems_web_low[, 1:npems_lat, drop = F]
        }
        if (npems_lat > ncol(pems_web_high)) {
          warning(paste("For web no.", i, ", only",  ncol(pems_web_high) , "PEMs will be used for the higher trophic level latent trait,\n since this network only has", ncol(webs[[i]]$web), "higher trophic level species."))
        } else {
          pems_web_high <- pems_web_high[, 1:npems_lat, drop = F]
        }
      }
    }
    webs[[i]]$pems <- list(low = pems_web_low, high = pems_web_high)
    webs[[i]]$abuns <- list(low = abun_low[[i]], high = abun_high[[i]])
  }
  if (is.null(traits_low)) {
    for (i in 1:length(networks)) webs[[i]]$traits <- NULL
  } else {
    for (i in 1:length(networks)) {
      traits_web_low <- traits_low[which(rownames(traits_low) %in% rownames(networks[[i]])), , drop = F]
      traits_web_low <- traits_web_low[order(rownames(traits_web_low)), , drop = F] # sort species alphabetically
      traits_web_high <- traits_high[which(rownames(traits_high) %in% colnames(networks[[i]])), , drop = F]
      traits_web_high <- traits_web_high[order(rownames(traits_web_high)), , drop = F] # sort species alphabetically
      webs[[i]]$traits <- list(low = traits_web_low, high = traits_web_high)
    }
  }
  
  out <- list(trees = trees, traits_all = traits_all, networks = webs)
  class(out) <- "tapnet"
  return(out)
}