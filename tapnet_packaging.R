library(devtools)
setwd("~/Data/aktuell/Networks/tapnet/tapnet")
document() # process R-functions into .RD files, change namespace
#tools::compactPDF("vignettes/", gs_quality = "ebook")

## Starting from here, the performance seems to vary:
setwd("..")
build("tapnet", args="--compact-vignettes=gs+qpdf", binary=F) 
install("tapnet") # install on the computer
devtools::check("tapnet", args="--as-cran")
# I don't get why it returns a warning about compactable PDF: it IS compacted and there is nothing anyone can do about this in the build!

## Using winbuilder (https://win-builder.r-project.org/upload.aspx) to check before submission, I repeatedly get funny problems when using the devtools-developed package. I suspect that some funny temporary file is involved, as the local errors point to var/local/cc/temp...
## Command line building works (mostly) fine.
R CMD build tapnet --compact-vignettes=gs+qpdf
R CMD check tapnet_0.4.tar.gz --as-cran

# %%%%%%%%%%%%%%%%%%%%%%%%%     %%%%%%%%%%%%%%%%%%%%%%%%%      %%%%%%%%%%%%%%%%%%%%%%%%%     %%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%     %%%%%%%%%%%%%%%%%%%%%%%%%      %%%%%%%%%%%%%%%%%%%%%%%%%     %%%%%%%%%%%%%%%%%%%%%%%%%
# fix things:
# use NON-STANDARDISED trait-matching function, so that Pmax is always 1, no matter which value sigma takes; prevents function to have a maximum at sigma=0
# re-run Tinoco-analysis: can delta be estimated properly (i.e. not as random value, but as optimum)? Do different optimisers yield same value?

# %%%%%%%%%%%%%%%%%%%%%%%%%     %%%%%%%%%%%%%%%%%%%%%%%%%      %%%%%%%%%%%%%%%%%%%%%%%%%     %%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%     %%%%%%%%%%%%%%%%%%%%%%%%%      %%%%%%%%%%%%%%%%%%%%%%%%%     %%%%%%%%%%%%%%%%%%%%%%%%%
# prepare things:
#save(humm_traits, humm_tree, networks, plant_traits, plant_tree, file="data/Tinoco.rda")

# test things:

library(tapnet)
data(Tinoco)
tap <- make_tapnet(tree_low = plant_tree, tree_high = humm_tree, networks = networks[2:3], traits_low = plant_traits, traits_high = humm_traits, abun_low=plant_abun[2:3], abun_high=humm_abun[2:3], npems_lat = 4)
fit <- fit_tapnet(tap, fit.delta = T) # fits on network 2 and 3
str(gof_tapnet(fit) )
pred1 <- predict_tapnet(fit, abuns=tap$networks[[1]]$abuns) # predict to forest network
cor(as.vector(pred1*sum(tap$networks[[1]]$web)), as.vector(tap$networks[[1]]$web)) # correlation with observed interactions
# slightly less bad with exp-link:
#cor(as.vector(exp(pred1)/sum(exp(pred1))*sum(tap$networks[[1]]$web)), as.vector(tap$networks[[1]]$web)) # correlation with observed interactions
fit

# test TmatchMatrixList option:
data(Tinoco)
tap2 <- make_tapnet(tree_low = plant_tree, tree_high = humm_tree, networks = networks[2], traits_low = plant_traits, traits_high = humm_traits, abun_low=plant_abun[2], abun_high=humm_abun[2], npems_lat = 4)
# make a mask:
Mask <- list(networks[[2]])
Mask[[1]][,] <- 1
Mask[[1]][10:11, 10:11] <- 0# 1e-8 # give unlikely events a VERY small value (less than 1/prod(dim(web)) ) 
fit2 <- fit_tapnet(tap2, TmatchMatrixList=Mask)
fit2 # different to without the mask




data(Tinoco)
tap <- make_tapnet(tree_low = plant_tree, tree_high = humm_tree, networks = networks[2:3], 
     traits_low = plant_traits, traits_high = humm_traits, npems_lat = 4)
system.time(fit <- fit_tapnet(tap, hessian=T, method="Nelder")) # fits to networks 2 and 3 only # with SANN: 252 s; with Nelder:  s
str(fit) 
(diag(solve(fit$opt$hessian)))

# test bjornloglik:
set.seed(2)
M <- matrix(rpois(80*100, 1), 80, 100)
P <- M/sum(M)
l <- sample.int(ncol(M)) # aims to test whether permutating entries makes difference to loglik (no!)
bjornloglik(M[,l], P[,l])

# test loglik_tapnet
fit_tapnet(tap, fit.delta=F, hessian=T) 
fit_tapnet(tap, fit.delta=F, obj_function="least squares", hessian=T) 
fit_tapnet(tap, fit.delta=F, obj_function="bjorn", hessian=T)
# all three yield rather different estimates; hessians look fine (for fitted parameters)
fit_tapnet(tap, fit.delta=F, obj_function="bjorn", hessian=T, tmatch_type_pem = "otto")  # should step with specific error message: works!

# test simulate_tapnet
si <- simulate_tapnet(10, 15, ntraits_nopem = 0, ntraits_pem = 0)
A <- si$networks[[1]]$abuns$low %*% t(si$networks[[1]]$abuns$high)
A <- A/sum(A)
cor(as.vector(A), as.vector(si$networks[[1]]$I_mat)) # this should be 1, as we use neither traits nor PEMs
# shows that simulate_tapnet uses PEMs even when ntraits_pem=0!!

#### wishlist: ####
* add function similar to select_relevant_pems to cut phylotree to species present EVEN when calling "use.all.pems=T" in make_tapnet (wish by Amanda)
* add correlation method option to gof (not only Spearman)
* add deOptim as optimiser (more robust but still fast, I hope)
* find data suitable for analysis and include in the package, e.g. https://doi.org/10.5061/dryad.nk98sf7sc (Wang et al., 2020). 
* add check for non-integer response?
* make a function for trait matching for phenology: give distribution for tmatch, e.g. Gaussian or uniform; or directly provide a phenology-match-matrix as input





#### use ranger for example ####
library(tapnet)
data(Tinoco)
# Note: use npems_lat = 105, use.all.pems=T to make sure ALL PEMs are being used; otherwise tapnet2df fails.
tap <- make_tapnet(tree_low = plant_tree, tree_high = humm_tree, networks = networks[2:3], traits_low = plant_traits, traits_high = humm_traits, abun_low=plant_abun[2:3], abun_high=humm_abun[2:3], npems_lat = 105, use.all.pems=T)
tapdats <- tapnet2df(tapnetObject=tap)
head(tapdats, 3)
tap1 <- make_tapnet(tree_low = plant_tree, tree_high = humm_tree, networks = networks[1], traits_low = plant_traits, traits_high = humm_traits, abun_low=plant_abun[1], abun_high=humm_abun[1], npems_lat = 105, use.all.pems=T)
tapdats1 <- tapnet2df(tapnetObject=tap1)
# ranger
library(ranger)
frf <- ranger(interactions ~., data=tapdats[, -c(1,2,4)])
plot(predict(frf, data=tapdats1)$predictions+1, tapdats1$interactions+1, log="xy")
abline(0,1)
cor(predict(frf, data=tapdats1)$predictions+1, tapdats1$interactions+1)

# mlp
library(RSNNS)
fmlp <- mlp(tapdats[, -c(1:4)], tapdats[, 3], size=c(10, 10, 10, 10))
predict(fmlp, newdata=tapdats1[, -c(1:3)])
#hm; try cito for better defaults?



#### modifying and trialling TAPNET ####

fit_tapnet -> optim -> logLik -> simnetfromtap


#### run without PEMs ####
library(tapnet)
data(Tinoco)
tap <- make_tapnet(tree_low = plant_tree, tree_high = humm_tree, networks = networks[3], traits_low = plant_traits, traits_high = humm_traits, abun_low=plant_abun[3], abun_high=humm_abun[3], npems_lat = 1)
# fit <- fit_tapnet(tap, fit.delta=F) # fits on network 1, only if npems_lat > 0!!
#fit <- fit_tapnet(tap, fit.delta=F, lambda=1) # fits on network 1
source("tapnet/R/helper_tapnet.R")
source("tapnet/R/simnetfromtap.R")
library(MPSEM)
source("tapnet/R/make_tapnet.R")
source("tapnet/R/fit_tapnet.R")
fit <- fit_tapnet(tap, fit.delta=F, tmatch_type_pem = "no") # fits on network 1, does not use PEMs
fit


tapnet=tap; fit.delta=F; lambda=0; tmatch_type_pem="normal"; ini = NULL; tmatch_type_obs = "normal"# for fit_tapnet
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

  
  
#### trial embedding ####

library(tapnet)
data(Tinoco)
# Note: use npems_lat = 105, use.all.pems=T to make sure ALL PEMs are being used; otherwise tapnet2df fails.
tap <- make_tapnet(tree_low = plant_tree, tree_high = humm_tree, networks = networks[2:3], traits_low = plant_traits, traits_high = humm_traits, abun_low=plant_abun[2:3], abun_high=humm_abun[2:3], npems_lat = 105, use.all.pems=T)
tapdats <- tapnet2df(tapnetObject=tap)
head(tapdats, 3)
tap1 <- make_tapnet(tree_low = plant_tree, tree_high = humm_tree, networks = networks[1], traits_low = plant_traits, traits_high = humm_traits, abun_low=plant_abun[1], abun_high=humm_abun[1], npems_lat = 105, use.all.pems=T)
tapdats1 <- tapnet2df(tapnetObject=tap1)

# using t-stochastic neighbourhood embedding:
library(Rtsne)
emb <- Rtsne(tapdats[,-c(1,2,4)])
plot(emb$Y, asp=1)
# hm; no way to predict from this
# setting input to NA leads to error

# using node2vec: network information only, no covariates!
library(node2vec)
# make edge-with-weights object, removing all 0 links:
x <- tapdats[which(tapdats[,3] > 0), c(1,2,3)]
unique(x[,2]) # 4 letters needed to uniquely identify species
x[, 1] <- sapply(x[,1], substr, 1, 4)
x[, 2] <- sapply(x[,2], substr, 1, 4)
x[, 3] <- log(x[,3] + 1)
head(x)

emb <- node2vecR(x, walk_length=10, p=2, dim=5)
emb



#### trial neural networks ####

library(tapnet)
data(Tinoco)
# Note: use npems_lat = 105, use.all.pems=T to make sure ALL PEMs are being used; otherwise tapnet2df fails.
tap <- make_tapnet(tree_low = plant_tree, tree_high = humm_tree, networks = networks[2:3], traits_low = plant_traits, traits_high = humm_traits, abun_low=plant_abun[2:3], abun_high=humm_abun[2:3], npems_lat = 105, use.all.pems=T)
tapdats <- tapnet2df(tapnetObject=tap)
head(tapdats, 3)
tap1 <- make_tapnet(tree_low = plant_tree, tree_high = humm_tree, networks = networks[1], traits_low = plant_traits, traits_high = humm_traits, abun_low=plant_abun[1], abun_high=humm_abun[1], npems_lat = 105, use.all.pems=T)
tapdats1 <- tapnet2df(tapnetObject=tap1)


library(nnet)
fnn <- nnet(x=tapdats[, -c(1:5)], y=tapdats$interactions, linout=T, decay=0.00001, size=20)
preds <- predict(fnn, newdata=tapdats1[, -c(1:5)])
plot(preds, tapdats1$interactions)
abline(0,1)


library(neuralnet)
fnet <- neuralnet(interactions ~ . , data=tapdats[, -c(1,2,4)], hidden=c(10, 7, 10), linear.output=T, stepmax=1e6)
str(fnet)
predict(fnet, newdata=tapdats1[, -c(1,2,4)])

 


#### trait matching as in Pichler et al. (2019) ####
# https://github.com/MaximilianPi/TraitMatching
library(tapnet)
data(Tinoco)
# Note: use npems_lat = 105, use.all.pems=T to make sure ALL PEMs are being used; otherwise tapnet2df fails.
tap <- make_tapnet(tree_low = plant_tree, tree_high = humm_tree, networks = networks[2:3], traits_low = plant_traits, traits_high = humm_traits, abun_low=plant_abun[2:3], abun_high=humm_abun[2:3], npems_lat = 105, use.all.pems=T)

devtools::install_github(repo = "https://github.com/MaximilianPi/TraitMatching", subdir = "TraitMatching")
library(TraitMatching)
#sim = simulateInteraction(weights = list(main = 0.1, inter = c(0.3,0.3,0.3)))
# sim$A, sim$B, sim$binar()
#community = createCommunity(A, B, Z)
# A and B are species groups with traits, respectively; Z is observed interactions
A <- cbind(abun=tap$networks[[1]]$abuns[[1]], tap$networks[[1]]$traits[[1]], tap$networks[[1]]$pems[[1]][names(tap$networks[[1]]$abuns[[1]]),])
B <- cbind(abun=tap$networks[[1]]$abuns[[2]], tap$networks[[1]]$traits[[2]], tap$networks[[1]]$pems[[2]][names(tap$networks[[1]]$abuns[[2]]),])
Z <- tap$networks[[1]]$web
comm <- createCommunity(a=A, b=B, z=Z, impute=F, log=T)
# does not work; creates empty list and I cannot be bothered to fix it
#fit = runTM(community = comm, method = "RF", iters = 20L)
