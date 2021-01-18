#' Simulates a network from parameters, abundances, traits and phylogeny provided
#'
#' The workhorse function of this package, called by various other functions to construct a bipartite interaction network
#'
#' Details in here!
#'
#' @aliases simnetfromtap
#' 
#' @param traits a named ("low"/"high") list of species-named trait data matrices for lower and higher trophic level;
#' @param abuns  a named ("low"/"high") list of species-named abundance vectors for lower and higher trophic level;
#' @param paramsList a list of parameter values with six elements: [[1]] and [[2]]: two vectors of linear combination parameters (importance values, one vector for each trophic level); [[3]]: a single shift parameter added to linear combination of higher trophic level PEMs; [[4]]: a single trait matching parameter for the PEMs; [[5]]: a vector of trait matching parameters for observed traits; [[6]]: a non-negative scalar delta to weight the importance of abundances 
#' @param pems  a named ("low"/"high") list of two species-named data frames (PEMs of lower and higher trophic level);
#' @param tmatch_type_pem type of trait matching function for latent traits (any name accepted by \code{\link{tmatch}}, currently "normal" and "shiftlnorm");
#' @param tmatch_type_obs  type of trait matching function for observed traits (see previous argument).
#' 
#' @return A named interaction matrix.
#' 
#' @references Benadi et al. in prep
#'
#' @author Gita Benadi <gita.benadi@biom.uni-freiburg.de> and Carsten Dormann <carsten.dormann@biom.uni-freiburg.de>
#'
#' @examples
#' \dontrun{
#' # whatever is in here
#' }
#'
#' @export
simnetfromtap <- function(traits, # named list of trait data matrices for lower and higher trophic level
                          abuns, # named list of abundance vectors for lower and higher trophic level
                          paramsList, # list of parameter values with six elements:
                          # [[1]] and [[2]]: two vectors of linear combination parameters
                          # (importance values, one vector for each trophic level)
                          # [[3]]: a single shift parameter added to linear combination of higher trophic level PEMs
                          # [[4]]: a single trait matching parameter for the PEMs
                          # [[5]]: a vector of trait matching parameters for observed traits
                          # [[6]]: a non-negative scalar delta to weight the importance of abundances 
                          pems, #  a list of two data frames (PEMs of lower and higher trophic level)
                          tmatch_type_pem, # type of trait matching function for latent traits
                          tmatch_type_obs  # type of trait matching function for observed traits
) {
  # Simulate interaction networks based on traits, abundances and phylogenies ("tap")
  
  # 0. alphabetically sort all entries!!
  if (!is.null(traits$low)) traits$low <- traits$low[order(rownames(traits$low)), , drop=FALSE]
  if (!is.null(traits$high)) traits$high <- traits$high[order(rownames(traits$high)), , drop=FALSE]
  
  abuns$low <- abuns$low[order(names(abuns$low))]
  abuns$high <- abuns$high[order(names(abuns$high))]
  
  pems$low <- pems$low[order(rownames(pems$low)), , drop=FALSE]
  pems$high <- pems$high[order(rownames(pems$high)), , drop=FALSE]
  
  # 1. Matching of observed traits
  if (is.null(traits$high) | is.null(traits$low)){
    # if we simulate a network without observed traits, we set T to a matrix full of 1s:
    T_mat <- matrix(1, nrow = length(abuns$low), ncol = length(abuns$high))
    T_mat <- T_mat / sum(T_mat)
  } else {
    T_mat <- tmatch(t(outer(traits$high[, 1], traits$low[, 1], "-")), type = tmatch_type_obs[1],
                    width = paramsList[[5]][1])
    rownames(T_mat) <- rownames(traits$low)
    colnames(T_mat) <- rownames(traits$high)
    T_mat <- T_mat / sum(T_mat)
    if (ncol(traits$low) > 1) {
      if (length(tmatch_type_obs) == 1) tmatch_type_obs <- rep(tmatch_type_obs, ncol(traits$low))
      if (length(paramsList[[5]]) == 1) paramsList[[5]] <- rep(paramsList[[5]], ncol(traits$low))
      for (i in 2:ncol(traits$low)) {
        T_mat_next <- tmatch(t(outer(traits$high[, i], traits$low[, i], "-")),
                             type = tmatch_type_obs[i], width = paramsList[[5]][i])
        rownames(T_mat) <- rownames(traits$low)
        colnames(T_mat) <- rownames(traits$high)
        T_mat_next <- T_mat_next / sum(T_mat_next)
        T_mat <- T_mat * T_mat_next
      }
    }
  }
  
  # 2. Compute linear combinations of PEMs:
  nspec_low <- length(abuns$low)
  nspec_high <- length(abuns$high)
  a_mat_low <- matrix(rep(paramsList[[1]], nspec_low), nrow = nspec_low, byrow = TRUE)
  a_mat_high <- matrix(rep(paramsList[[2]], nspec_high), nrow = nspec_high, byrow = TRUE)
  lat_low <- as.vector(scale(rowSums(a_mat_low * pems[[1]])))
  lat_high <- as.vector(scale(rowSums(a_mat_high * pems[[2]]))) + paramsList[[3]] # Apply shift to higher trophic level
  
  # Compute interaction probabilities based on linear combinations:
  L_mat <- tmatch(t(outer(lat_high, lat_low, "-")), type = tmatch_type_pem, width = paramsList[[4]])
  rownames(L_mat) <- rownames(pems[[1]])
  colnames(L_mat) <- rownames(pems[[2]])
  L_mat <- L_mat / sum(L_mat) # scale to sum of 1
  
  # Interaction probabilities based on abundances:
  A_mat <- as.matrix(abuns$low) %*% t(abuns$high)
  A_mat <- A_mat / sum(A_mat) # scale to sum of 1
  
  # sort A_mat:
  #	A_mat <- A_mat[order(rownames(A_mat)), order(colnames(A_mat))]
  
  # T_mat <- T_mat / sum(T_mat)
  # note that T_mat is not scaled, allowing the optimiser to adjust balance between TRUE and L!
  
  if (is.null(paramsList[["delta"]])) {
    LT <- (L_mat * T_mat) / sum(L_mat * T_mat)
  } else {
    delta <- paramsList[["delta"]]
    #A_mat <- (A_mat^plogis(delta)) / sum(A_mat^plogis(delta)) # renormalise to sum to 1
    #L_mat <- (L_mat^plogis(delta)) / sum(L_mat^plogis(delta)) # renormalise to sum to 1
    LT <- (L_mat * T_mat)^(plogis(delta)) / sum((L_mat * T_mat)^(plogis(delta))) # delta now balances abundance vs. traits
  }
  
  # Put everything together:
  I_mat <- A_mat * LT #L_mat * T_mat # c allows adjustment of importance of abundances
  I_mat <- I_mat / sum(I_mat)  # scaling the matrix to a sum of 1
  
  return(I_mat)
}
