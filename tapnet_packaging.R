library(devtools)
setwd("~/Data/aktuell/Networks/tapnet/tapnet")
document() # process R-functions into .RD files, change namespace
#tools::compactPDF("vignettes/", gs_quality = "ebook")

## Starting from here, the performance seems to vary:
setwd("..")
build("tapnet", args="--compact-vignettes=gs+qpdf", binary=F) 
#install("tapnet") # install on the computer
devtools::check("tapnet", args="--as-cran")
# I don't get why it returns a warning about compactable PDF: it IS compacted and there is nothing anyone can do about this in the build!

## Using winbuilder (https://win-builder.r-project.org/upload.aspx) to check before submission, I repeatedly get funny problems when using the devtools-developed package. I suspect that some funny temporary file is involved, as the local errors point to var/local/cc/temp...
## Command line building works (mostly) fine.
R CMD build tapnet --compact-vignettes=gs+qpdf
R CMD check tapnet_0.4.tar.gz --as-cran

# %%%%%%%%%%%%%%%%%%%%%%%%%     %%%%%%%%%%%%%%%%%%%%%%%%%      %%%%%%%%%%%%%%%%%%%%%%%%%     %%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%     %%%%%%%%%%%%%%%%%%%%%%%%%      %%%%%%%%%%%%%%%%%%%%%%%%%     %%%%%%%%%%%%%%%%%%%%%%%%%
# prepare things:
#save(humm_traits, humm_tree, networks, plant_traits, plant_tree, file="data/Tinoco.rda")

# test things:

library(tapnet)
data(Tinoco)
tap <- make_tapnet(tree_low = plant_tree, tree_high = humm_tree, networks = networks[2:3], traits_low = plant_traits, traits_high = humm_traits, abun_low=plant_abun[2:3], abun_high=humm_abun[2:3], npems_lat = 4)
fit <- fit_tapnet(tap, fit.delta = T) # fits on network 2 and 3
gof_tapnet(fit) 
pred1 <- predict_tapnet(fit, abuns=tap$networks[[1]]$abuns) # predict to forest network
cor(as.vector(pred1*sum(tap$networks[[1]]$web)), as.vector(tap$networks[[1]]$web)) # correlation with observed interactions
fit



data(Tinoco)
tap <- make_tapnet(tree_low = plant_tree, tree_high = humm_tree, networks = networks[2:3], 
     traits_low = plant_traits, traits_high = humm_traits, npems_lat = 4)
system.time(fit <- fit_tapnet(tap, hessian=T, method="Nelder")) # fits to networks 2 and 3 only # with SANN: 252 s; with Nelder:  s
str(fit) 
(diag(solve(fit$opt$hessian)))


# wishlist:
* add correlation method option to gof (not only Spearman)
* add deOptim as optimiser (more robust but still fast, I hope)
* find data suitable for analysis and include in the package, e.g. https://doi.org/10.5061/dryad.nk98sf7sc (Wang et al., 2020). 




#### trialling TAPNET ####

fit_tapnet -> optim -> logLik -> simnetfromtap


#### run without PEMs ####
library(tapnet)
data(Tinoco)
tap <- make_tapnet(tree_low = plant_tree, tree_high = humm_tree, networks = networks[3], traits_low = plant_traits, traits_high = humm_traits, abun_low=plant_abun[3], abun_high=humm_abun[3], npems_lat = 0)
# fit <- fit_tapnet(tap, fit.delta=F) # fits on network 1, only if npems_lat > 0!!
#fit <- fit_tapnet(tap, fit.delta=F, lambda=1) # fits on network 1
source("tapnet/R/helper_tapnet.R")
source("tapnet/R/simnetfromtap.R")
library(MPSEM)
source("tapnet/R/make_tapnet.R")
source("tapnet/R/fit_tapnet.R")
fit <- fit_tapnet(tap, fit.delta=F, tmatch_type_pem = "no") # fits on network 1, does not use PEMs
fit


tapnet=tap; fit.delta=F; lambda=0; tmatch_type_pem="no"; ini = NULL; tmatch_type_obs = "normal"# for fit_tapnet
# now go through fit_tapnet up to line 92 (before optim) to produce all variables for use in simnetfromtap
networks <- tap$networks
# now everything is set up to run the loglik for different values of sigma (see below)!

#params = ini; method = "Nelder"; maxit=500; obj_function = "multinom"; hessian=T # now go through logLik, up to line 170 (which calls simnetfromtap), then switch to there
#networks <- tap$networks; traits = tap$networks[[i]]$traits; abuns = tap$networks[[i]]$abuns; pems = tap$networks[[i]]$pems
#tmatch_type_pem = tmatch_type_pem, tmatch_type_obs = tmatch_type_obs




# evaluate logLik for different values of tmatch_width_obs1:
trials <- ini

ell=0
j=1
sigma <- seq(0.1, 10, len=100)
for (i in sigma){
  trials[2] <- i
  ell[j] <- loglik_tapnet(params=trials, networks=networks, tmatch_type_pem = "no", tmatch_type_obs = "normal", fit.delta=F)
  j = j + 1
}
plot(sigma, ell, type="l") # fine for Tinoco 1 and 2
plot(sigma, ell, type="l", ylim=c(2000, 4000)) # very weak optimum for network 3 !!
abline(v=2.77)



#### run fit_tapnet with different random starting values ####



3. make matrix with high matching: does it show a clear optimum when evaluating the loglik?
