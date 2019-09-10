#' Convert tapnet object into data.frame
#'
#' Function allows direct use of data prepared for tapnet analysis by other statistical methods, e.g. regression approaches
#'
#' This function simply puts all data into a data.frame, with each row an entry in the network matrix.
#'
#' @aliases tapnet2df
#' 
#' @param fit results of applying fit_tapnet to the tapnet object;
#' 
#' @return A data.frame containing network observations, PEMs, traits and abundances for regression-type analysis. 
#' 
#' @references Benadi et al. in prep
#'
#' @author Carsten Dormann <carsten.dormann@biom.uni-freiburg.de>
#'
#' @examples
#  ex <- simulate_tapnet(nlower=10, nhigher=50, ntraits_pem=3, ntraits_nopem=2, Nwebs = 3)
#' df <- tapnet2df(ex)
#' head(df)
#' \dontrun{
#'   library(ranger)
#'   frf <- ranger(interactions ~ ., data=df[, -c(1:2)], importance="impurity")
#'   sort(importance(frf), decreasing=TRUE)
#' }
#'
#' @export
tapnet2df <- function(tapnetObject){
  # function to convert a tapnet-object into a data.frame to be used, e.g., by any typical machine-learning approach
  
  # TODO: fill all networks with species so that all networks contain all species; that may make the use easier for some other functions
  # idea: duplicate tapnetObject and replace existing networks with full networks, trait matrices and PEMs
  
  
  if (!inherits(tapnetObject, "tapnet")) stop("This is no tapnet object! Please use function 'make_tapnet' or 'simulate_tapnet' to create it.")
  
  # put a large for-loop around the following code (run through networks)
  #extract elements: 
  web <- tapnetObject$networks
  # if (length(web) > 1){
  #	warning("Currently tapnet2df can only handle a single network. The first in the list of networks provided will be used.")
  #	one.network <- tapnetObject$networks[[1]]
  #} 
  
  for (i in seq_along(web)){ # extract information for each network
    
    one.network <- tapnetObject$networks[[i]]
    
    # stack interactions:
    dats <- stack(as.data.frame(one.network$web))
    colnames(dats) <- c("interactions", "IDhigher")
    dats$IDlower <- rep(rownames(one.network$web), ncol(one.network$web))
    if (length(web) > 1) dats$webNumber <- i
    
    # extract PEMs from tapnet object and add pems for lower level:
    blubb <- one.network$pems$low
    colnames(blubb) <- paste0("pemL", colnames(blubb))
    blubb$names <- rownames(blubb)
    dats <- merge(dats, blubb, by.x="IDlower", by.y="names")
    rm(blubb)
    
    # extract PEMs from tapnet object and add pems for higher level:
    blubb <- one.network$pems$high
    colnames(blubb) <- paste0("pemH", colnames(blubb))
    blubb$names <- rownames(blubb)
    dats <- merge(dats, blubb, by.x="IDhigher", by.y="names")
    rm(blubb)
    
    # add abundances lower
    if (is.null(one.network$abun$low)){
      warning("No lower abundances provided, marginal totals will be used.")
      lowerAbun <- data.frame("abunL" = rowSums(one.network$web))
      lowerAbun$names <- rownames(lowerAbun)
    } else { # take values from argument:
      lowerAbun  <- data.frame("abunL" = one.network$abun$low)
      lowerAbun$names <- rownames(one.network$web)
    }
    dats <- merge(dats, lowerAbun, by.x = "IDlower", by.y = "names") # merge in lower level IDs
    rm(lowerAbun)
    
    # add abundances higher
    if (is.null(one.network$abun$high)){
      warning("No higher abundances provided, marginal totals will be used.")
      higherAbun <- data.frame("abunH" = colSums(one.network$web))
      higherAbun$names <- rownames(higherAbun)
    } else { # take values from argument:
      higherAbun  <- data.frame("abunH" = one.network$abun$high)
      higherAbun$names <- colnames(one.network$web)
    }
    dats <- merge(dats, higherAbun, by.x = "IDhigher", by.y = "names") # merge in lower level IDs
    rm(higherAbun)
    
    # add observed traits lower:
    if (is.null(one.network$traits$low)){
      warning("No values for lower-level observed traits provided and hence none used.")
    } else {
      # add traits for lower level:
      lowerTraits <- as.data.frame(tapnetObject$traits$low)
      colnames(lowerTraits) <- paste0("traitL", colnames(lowerTraits))
      lowerTraits$names <- rownames(lowerTraits)#rownames(one.network$web)
    }
    dats <- merge(dats, lowerTraits, by.x = "IDlower", by.y = "names")
    rm(lowerTraits)
    
    # add observed traits higher:
    if (is.null(one.network$traits$high)){
      warning("No values for higher-level observed traits provided and hence none used.")
    } else {
      # add traits for higher level:
      higherTraits <- as.data.frame(tapnetObject$traits$high)
      colnames(higherTraits) <- paste0("traitH", colnames(higherTraits))
      higherTraits$names <- rownames(higherTraits)# colnames(one.network$web)
    }
    dats <- merge(dats, higherTraits, by.x = "IDhigher", by.y = "names")
    rm(higherTraits)
    rm(one.network) 
    
    # store each network underneath the previous network:
    if (i == 1) allWebs <- dats else allWebs <- rbind(allWebs, dats)
    # TODO: this actually should be done by a merge, as each network may have different number of PEMS and traits!
    # so: transpose, merge, re-transpose (I guess)
    
    rm(dats)
  } # end for loop through webs  
  rm(web)
  
  return(allWebs)
} 
