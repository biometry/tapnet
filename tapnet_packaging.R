library(devtools)
setwd("~/Data/aktuell/Networks/tapnet/tapnet")
document() # process R-functions into .RD files, change namespace
setwd("..")
build("tapnet", args="--compact-vignettes")   # build the .tar.gz file  ## takes 11 min because of the processing of the .Rnw file!!
#build("tapnet", args="--no-build-vignettes", clean_doc=FALSE) # if vignette hasn't changed; does not work for package to be uploaded to CRAN!!
install("tapnet") # install on the computer
# check_built("tapnet") # is included in:
devtools::check("tapnet")

# hm; I get an error message from winbuilder with this, saying: 
* checking index information ... WARNING
Vignettes with missing or empty \VignetteIndexEntry:
  CaseStudy.Rnw
See sections 'The INDEX file' and 'Package subdirectories' in the
'Writing R Extensions' manual.
# And indeed, the installed tapnet package html has no link to the CaseStudy. 
# When using the command line, the link is there ...
R CMD build tapnet --compact-vignettes
R CMD check tapnet_0.3.tar.gz --as-cran




# prepare things:
#save(humm_traits, humm_tree, networks, plant_traits, plant_tree, file="data/Tinoco.rda")

# test things:

library(tapnet)
data(Tinoco)
tap <- make_tapnet(tree_low = plant_tree, tree_high = humm_tree, networks = networks[2:3], traits_low = plant_traits, traits_high = humm_traits, abun_low=plant_abun[2:3], abun_high=humm_abun[2:3], npems_lat = 4)
fit <- fit_tapnet(tap) # fits on network 2 and 3
gof_tapnet(fit) 
pred1 <- predict_tapnet(fit, abuns=tap$networks[[1]]$abuns) # predict to forest network
cor(as.vector(pred1*sum(tap$networks[[1]]$web)), as.vector(tap$networks[[1]]$web)) # correlation with observed interactions





data(Tinoco)
tap <- make_tapnet(tree_low = plant_tree, tree_high = humm_tree, networks = networks[2:3], 
     traits_low = plant_traits, traits_high = humm_traits, npems_lat = 4)
system.time(fit <- fit_tapnet(tap, hessian=T, method="Nelder")) # fits to networks 2 and 3 only # with SANN: 252 s; with Nelder:  s
str(fit) 
(diag(solve(fit$opt$hessian)))


# wishlist:
* add correlation method option to gof (not only Spearman)