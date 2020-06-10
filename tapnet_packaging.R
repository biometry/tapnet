library(devtools)
setwd("~/Data/aktuell/Networks/tapnet/")
document() # process R-functions into .RD files, change namespace
setwd("..")
build("tapnet", path="tapnet")   # build the .tar.gz file
install("tapnet") # install on the computer

# check_built("tapnet") # is included in:^
devtools::check("tapnet")


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
