
## testtest: start with yesterdays tapnet, remove vignette, re-build and re-submit to winbuilder:
fine!
## changed description: no author, no maintainer, suggest only knitr, vignetteBuilder knit and utils
fine!
## Re-run roxygen ("document") on vignette-free package: does it screw it up?
fine!
## Add vignette folder, roxygenise
!! error message about LOCK-file in deleted temp-folder (private/var/folders/...); Re-installed R, no effect
## copied all old R and Description function from yesterday onto current folder; re-submit to winbuilder
fine! (except for pdf compaction)
## added --compact-vignettes=gs+qpdf to tapnet-mainline, resubmit to release
fine! (except for PDF compaction)
## fixed error (testtest; deleted "build" folder), build package, submitted to develop! 
fine!
## mainline tapnet build with devtools and uploaded to release with two text snippets added to tapnet-package (e.g. reference)
fine!
## same mainline/devtools to devel:

  
# ideen: 
* utils in description?
* LF line endings corrupted (why?)
* roxygen screwed up

## winbuilder screws up: cache issues? --> new try on old-release with devtoolsbuild; no, same problem there
## removed vignette (in testtest) and uploaded to release: same problem
## release: packaged with devtools --> error:
** R
Error in parse(outFile) : 
  d:/temp/Rtmp0eBoln/R.INSTALL102386d2e5c71/tapnet/R/tapnet-package.R:1:1: unexpected '<'
1: <
  ^
  ERROR: unable to collate and parse R files for package 'tapnet'
## devel: packaged on ubuntu, command line: same error



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