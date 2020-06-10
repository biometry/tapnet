#' Hummingbird-flower networks
#' 
#' An example dataset for tapnet analysis
#'
#' These data are from the supplement of Tinoco et al. (2017) and contain three observed networks (forest at Mazan, shrubland at Llaviuco, cattle farm at Nero, all in Ecuador), along with traits of flowers and birds (corolla and beak length, respectively) as well as phylogenies and external abundances for all species. These data are in several ways special, but most of all because of the very high sampling effort that went into the networks.
#' 
#' For sake of clarity, we provide the data as separate objects. So when "Tinoco" is called, it will load seven objects: networks, humm_traits, humm_tree, humm_abun, plant_traits, plant_tree and plant_abun. To combine them into a useable tapnet object, use \code{\link{make_tapnet}}. Phylogenetic trees are of class "phylo" (as used/produced by \pkg{phytools}). Abundance data were provided independently of the other data directly by Boris (the other data are on dryad \url{http://dx.doi.org/10.5061/dryad.j860}). For the external abundances of hummingbirds, 12 point counts were performed in the same habitats where hummingbird - plant interactions were observed. "Abundance were obtained by averaging the abundance of each species per point count across the study period." "Plant abundances are averages across the study period." Many thanks to Boris for making his data freely available!
#'
#' @name Tinoco
#' @aliases humm_traits humm_tree  humm_abun networks plant_tree plant_traits plant_abun
#' @docType data
#' @references Tinoco, B. A.; Graham, C. H.; Aguilar, J. M. & Schleuning, M. Effects of hummingbird morphology on specialization in pollination networks vary with resource availability. \emph{Oikos} \bold{126}, 52-–60
#' 
#' @usage data(Tinoco)
#' @author Boris A. Tinoco Molina \email{btinoco@uazuay.edu.ec} collected the data; Carsten F. Dormann \email{carsten.dormann@biom.uni-freiburg.de} packaged them
#' @keywords data
"humm_traits"
"humm_trees"
