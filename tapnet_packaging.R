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
R CMD check tapnet_0.3.tar.gz --as-cran

# %%%%%%%%%%%%%%%%%%%%%%%%%     %%%%%%%%%%%%%%%%%%%%%%%%%      %%%%%%%%%%%%%%%%%%%%%%%%%     %%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%     %%%%%%%%%%%%%%%%%%%%%%%%%      %%%%%%%%%%%%%%%%%%%%%%%%%     %%%%%%%%%%%%%%%%%%%%%%%%%
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
* add deOptim as optimiser (more robust but still fast, I hope)
* find data suitable for analysis and include in the package, e.g. https://doi.org/10.5061/dryad.nk98sf7sc (Wang et al., 2020). 




#### old testing experiences .... ####
