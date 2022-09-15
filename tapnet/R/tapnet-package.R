#' tapnet: A package for fitting, and predicting from, bipartite interaction networks.
#'
#' The package provides three categories of important functions:
#' making tapnet objects, fitting the tapnet, and predicting to new data. All other functions are helpers or quality monitoring. 
#' The idea of using abundances, observed traits and phylogeny-derived latent traits to quantitatively fit and predict interaction 
#' is detailed in Benadi et al. (2021). This package accompanies the paper.
#' 
#' @section Making a tapnet object:
#' The typical tapnet object contains information on (1) an interaction network in the form of one or more matrices; (2) a phylogenetic tree for both groups; (3) observed traits for both groups, arranged in such a way that they are supposed to match when the difference is 0. These information can come separately, and the make_tapnet function aims at putting them into a tapnet object, as is expected by the functions in this package.
#' Also, the function simulate_tapnet can be used to simulate data of this format.
#' 
#' @section Fitting a tapnet:
#' Fitting of a tapnet refers to the idea of trying to find a function that predicts interactions depending on species abundances, their observed traits (and the match of these traits) and the match of unobserved "latent" traits constructed from the phylogeny. This fitted function can have quite a few parameter (we refer to the paper here), and fit_tapnet is the function to call all relevant optimisers and likelihood functions. 
#' At present, the observed interactions are evaluated as multinomial likelihood given the expectation from the parameters that are fitted. Also a least-square option is available. The fit can be evaluated in different ways, e.g. using the gof_tapnet function.
#' 
#' @section Predicting from a tapnet:
#' At present only predictions to altered abundances are possible. To predict to a new species, this means fitting the original network with this species in it, but no observed interactions. For prediction, the new abundances and the fitted tapnet need to be provided. See predict_tapnet for an example.
#' 
#' @section News/versions:
#' \describe{
#'   \item{0.4: 15-Sep-2022}{Added the option to fit networks without phylogenetic information. To do so, use \option{tmatch_type_pems="no"}. Rather experimental at this stage.}
#'   \item{0.3; 10-Jun-2019}{Added option "tapnet" to predict_tapnet in order to be able to use it on simulations. Cleaned up the code according to devtools::check-report.}
#'   \item{0.2; 19-Sep-2019}{Tinoco-data updated to now also contain external abundances (thanks Boris for providing these!). Help file for these data updated accordingly. By now we have used tapnet on many simulations and some real data and are fairly confident that the functions work as they should. The real test comes when data are not carefully prepared, with corrupted names and in non-alphabetical order, networks have NAs and so forth. Hopefully we have coded properly for these cases, too.}
#'   \item{0.1}{Initial version with all functions complete and a data set for demonstration. Code written by Gita Benadi, Carsten Dormann and Jochen Fründ, data provided by Boris Tinoco.}
#' }
#' 
#'@section References:
#' Benadi, G., Dormann, C.F., Fründ, J., Stephan, R.\& Vázquez, D.P. (2022) Quantitative prediction of interactions in bipartite networks based on traits, abundances, and phylogeny. \emph{The American Naturalist} \bold{199}, 841--854.
#'
#' @import stats
#' @importFrom ape keep.tip
#' @importFrom bipartite networklevel
#' @importFrom MPSEM Phylo2DirectedGraph PEM.build
#' @importFrom phytools pbtree
#' @importFrom utils stack
#' @importFrom vegan vegdist
#'
#' @docType package
#' @name tapnet-package
NULL