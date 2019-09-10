#' Hummingbird-flower networks
#' 
#' An example dataset for tapnet analysis
#'
#' These data are from the supplement of Tinoco et al. (2017) and contain three observed networks (forest, shrubland, cattle farm in Ecuador), along with traits of flowers and birds (corolla and beak length, respectively) as well as phylogenies for all species. These data are in several ways special, but most of all because of the very high sampling effort that went into the networks. Note that the data have no external abundances.
#' 
#' For sake of clarity, we provide the data as separate objects. So when "Tinoco" is called, it will load five objects: networks, humm_traits, humm_tree, plant_traits and plant_tree. To combine them into a useable tapnet object, use \code{\link{make_tapnet}}. Phylogenetic trees are of class "phylo" (as used/produced by \pkg{phytools}).
#'
#' @name Tinoco
#' @aliases humm_traits
#' @aliases humm_trees
#' @aliases networks
#' @aliases plant_trees
#' @aliases plant_traits
#' @docType data
#' @references Tinoco, B. A.; Graham, C. H.; Aguilar, J. M. & Schleuning, M. Effects of hummingbird morphology on specialization in pollination networks vary with resource availability. \emph{Oikos} \bold{126}, 52-–60
#' 
#' @usage data(Tinoco)
#' @author Carsten F. Dormann \email{carsten.dormann@biom.uni-freiburg.de}
#' @keywords data
"humm_traits"
"humm_trees"
