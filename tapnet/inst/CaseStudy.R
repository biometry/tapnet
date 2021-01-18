### R code from vignette source 'CaseStudy.Rnw'

###################################################
### code chunk number 1: setup
###################################################
library(knitr)
opts_chunk$set(fig.path='figure/', fig.align='center', fig.width=6, fig.height=6, fig.show='hold', cache=T, tidy=F, tidy.opts=list(width.cutoff=65), size="small")
#render_listings() # makes pretty formatting but allows lines to overflow
options(width = 80)


###################################################
### code chunk number 2: load data
###################################################
library(tapnet)
data(Tinoco)


###################################################
### code chunk number 3: construct tapnet objects
###################################################
# Produce tapnet objects using each network separately (1=Forest, 2=shrub, 3=cattle farm)
tapnet_web1 <- make_tapnet(tree_low = plant_tree, tree_high = humm_tree, 
                           networks = networks[1], traits_low = plant_traits, 
                           traits_high = humm_traits, abun_low=plant_abun[1], 
                           abun_high=humm_abun[1], npems_lat = 4)
tapnet_web2 <- make_tapnet(tree_low = plant_tree, tree_high = humm_tree, 
                           networks = networks[2], traits_low = plant_traits, 
                           traits_high = humm_traits, abun_low=plant_abun[2], 
                           abun_high=humm_abun[2], npems_lat = NULL)


###################################################
### code chunk number 4: CaseStudy.Rnw:129-135
###################################################
colnames(tapnet_web1$networks[[1]]$pems$low) # names of fitted PEMs
colnames(tapnet_web1$networks[[1]]$pems$high)

colnames(tapnet_web2$networks[[1]]$pems$low) # names of PEMs all present 
                                             
colnames(tapnet_web2$networks[[1]]$pems$high) # V_8 (high) is missing!


###################################################
### code chunk number 5: CaseStudy.Rnw:138-141
###################################################
tapnet_web2$networks[[1]]$pems$high$V_8 <- tapnet:::pems_from_tree(humm_tree)[colnames(
  tapnet_web2$networks[[1]]$web), "V_8"]
colnames(tapnet_web2$networks[[1]]$pems$high) # check: complete!


###################################################
### code chunk number 6: correlate PEMs and observed
###################################################
cor(cbind(tapnet_web1$networks[[1]]$pems$low, tapnet_web1$networks[[1]]$traits$low))
cor(cbind(tapnet_web1$networks[[1]]$pems$high, tapnet_web1$networks[[1]]$traits$high))


###################################################
### code chunk number 7: fit tapnet
###################################################
fit_web1 <- fit_tapnet(tapnet = tapnet_web1, method="SANN") # very slow, but reliable
#fit_web1 <- fit_tapnet(tapnet = tapnet_web1) # the default way
#fit_web1ln <- fit_tapnet(tapnet = tapnet_web1, tmatch_type_obs = "shiftlnorm", 
#                         ini=fit_web1$opt$par*2) # requires some tempering with ini
gof_web1_norm <- gof_tapnet(fit_web1)
gof_web1_norm


###################################################
### code chunk number 8: look at fit
###################################################
fit_web1


###################################################
### code chunk number 9: correlate observed and latent
###################################################
fitted_lin_low <- fit_web1$par_opt$lat_low[which(names(fit_web1$par_opt$lat_low) %in% 
                              colnames(tapnet_web1$networks[[1]]$pems$low))]
fitted_lat_low <- as.vector(scale(rowSums(matrix(fitted_lin_low, 
                            nrow = nrow(tapnet_web1$networks[[1]]$pems$low), 
                            ncol = ncol(tapnet_web1$networks[[1]]$pems$low), byrow = TRUE) * 
                            tapnet_web1$networks[[1]]$pems$low)))
cor(fitted_lat_low, tapnet_web1$networks[[1]]$traits$low)

fitted_lin_high <- fit_web1$par_opt$lat_high[which(names(fit_web1$par_opt$lat_high) %in% 
                                        colnames(tapnet_web1$networks[[1]]$pems$high))]
fitted_lat_high <- as.vector(scale(rowSums(matrix(fitted_lin_high, 
                              nrow = nrow(tapnet_web1$networks[[1]]$pems$high), 
                              ncol = ncol(tapnet_web1$networks[[1]]$pems$high), byrow = TRUE) * 
                              tapnet_web1$networks[[1]]$pems$high)))
cor(fitted_lat_high, tapnet_web1$networks[[1]]$traits$high)


###################################################
### code chunk number 10: correlate abundance and latent
###################################################
cor(fitted_lat_low, tapnet_web1$networks[[1]]$abuns$low)

cor(fitted_lat_high, tapnet_web1$networks[[1]]$abuns$high)


###################################################
### code chunk number 11: predict tapnet web2
###################################################
preds2.tapnet <- predict_tapnet(fit=fit_web1, abuns=tapnet_web2$networks[[1]]$abuns)
cor(as.vector(preds2.tapnet), as.vector(tapnet_web2$networks[[1]]$web))


###################################################
### code chunk number 12: predict tapnet web2 plot
###################################################
sum(tapnet_web2$networks[[1]]$web)
par(mar=c(5,5,1,1))
plot(preds2.tapnet*3979 + 1, tapnet_web2$networks[[1]]$web + 1, log="xy", las=1, 
  xlab="predicted number of interactions + 1", ylab="observed number of interactions + 1")
abline(0,1)


###################################################
### code chunk number 13: CaseStudy.Rnw:227-229
###################################################
dmultinom(as.vector(tapnet_web2$networks[[1]]$web), prob=as.vector(preds2.tapnet), 
          size=sum(tapnet_web2$networks[[1]]$web), log=T)


###################################################
### code chunk number 14: fit 2 networks example
###################################################
data(Tinoco)
tap <- make_tapnet(tree_low = plant_tree, tree_high = humm_tree, networks = networks[2:3], 
                   traits_low = plant_traits, traits_high = humm_traits, 
                   abun_low = plant_abun[2:3], abun_high=humm_abun[2:3] , npems_lat = 4)
fit <- fit_tapnet(tap) # uses two networks for fitting!
gof_tapnet(fit)
# predict to omitted forest network:
pred1 <- predict_tapnet(fit, abuns=list("low"=plant_abun[[1]], "high"=humm_abun[[1]] )) 
  
cor(as.vector(pred1*sum(networks[[1]])), as.vector(networks[[1]])) 


###################################################
### code chunk number 15: fit 2 networks example cont.
###################################################
tap <- make_tapnet(tree_low = plant_tree, tree_high = humm_tree, networks = networks[c(1,3)], traits_low = plant_traits, traits_high = humm_traits, abun_low = plant_abun[c(1,3)], 
                           abun_high=humm_abun[c(1,3)] , npems_lat = 4)
fit <- fit_tapnet(tap) # uses two networks for fitting!
pred1 <- predict_tapnet(fit, abuns=list("low"=plant_abun[[2]], "high"=humm_abun[[2]] )) # predict to omitted forest network
cor(as.vector(pred1*sum(networks[[2]])), as.vector(networks[[2]])) 


###################################################
### code chunk number 16: fit 2 networks example contcont.
###################################################
tap <- make_tapnet(tree_low = plant_tree, tree_high = humm_tree, networks = networks[1:2], traits_low = plant_traits, traits_high = humm_traits, abun_low = plant_abun[1:2], 
                           abun_high=humm_abun[1:2] , npems_lat = 4)
fit <- fit_tapnet(tap) # uses two networks for fitting!
pred1 <- predict_tapnet(fit, abuns=list("low"=plant_abun[[3]], "high"=humm_abun[[3]] )) # predict to omitted forest network
cor(as.vector(pred1*sum(networks[[3]])), as.vector(networks[[3]])) 


###################################################
### code chunk number 17: CaseStudy.Rnw:268-272
###################################################
preds2.abunonly <- (tapnet_web2$networks[[1]]$abuns$low / 
    sum(tapnet_web2$networks[[1]]$abuns$low)) %*% t(tapnet_web2$networks[[1]]$abuns$high / 
    sum(tapnet_web2$networks[[1]]$abuns$high)) * sum(tapnet_web2$networks[[1]]$web)
cor(as.vector(preds2.abunonly), as.vector(tapnet_web2$networks[[1]]$web))


###################################################
### code chunk number 18: abunonly plot
###################################################
par(mar=c(5,5,1,1))
plot(preds2.abunonly + 1, tapnet_web2$networks[[1]]$web + 1, log="xy", las=1, 
  xlab="predicted number of interactions + 1", ylab="observed number of interactions + 1")
abline(0,1)


###################################################
### code chunk number 19: CaseStudy.Rnw:281-283
###################################################
dmultinom(as.vector(tapnet_web2$networks[[1]]$web), prob=as.vector(preds2.abunonly), 
          size=sum(tapnet_web2$networks[[1]]$web), log=T)


###################################################
### code chunk number 20: make df from tapnet
###################################################
web1.df <- tapnet2df(tapnet_web1)
web2.df <- tapnet2df(tapnet_web2)
head(web1.df)


###################################################
### code chunk number 21: add traitmatch to df
###################################################
web1.df.extended <- cbind.data.frame(web1.df, "match"=(web1.df$traitHBill_length_mean_mm - 
                                                  web1.df$traitLCorolla_length_mm)^2 )
web2.df.extended <- cbind.data.frame(web2.df, "match"=(web2.df$traitHBill_length_mean_mm - 
                                                  web2.df$traitLCorolla_length_mm)^2 )


###################################################
### code chunk number 22: fit GAM
###################################################
library(mgcv)
gam2 <- gam(interactions ~ s(pemLV_1, pemHV_1, bs="ts", k=24) +s(pemLV_2, pemHV_3, bs="ts", 
           k=24) + s(traitLCorolla_length_mm, k=3) + s(traitHBill_length_mean_mm, k=3) + 
           s(match, k=3) + s(abunL, k=3) + s(abunH, k=3), data=web1.df.extended, family=nb, 
           gamma=1.4)
summary(gam2)
preds2.gam <- predict(gam2, newdata=web2.df.extended)
cor(exp(preds2.gam), web2.df$interactions)


###################################################
### code chunk number 23: predict 2 gam plot
###################################################
par(mar=c(5,5,1,1))
plot(exp(preds2.gam) +1 , web2.df$interactions + 1, log="xy", las=1, xlab="predicted 
     number of interactions + 1", ylab="observed number of interactions + 1")
abline(0,1)


###################################################
### code chunk number 24: CaseStudy.Rnw:333-335
###################################################
dmultinom(as.vector(tapnet_web2$networks[[1]]$web), prob=exp(preds2.gam) / 
            sum(exp(preds2.gam)), size=sum(tapnet_web2$networks[[1]]$web), log=T)


###################################################
### code chunk number 25: fit ranger web1
###################################################
library(ranger)
rf2 <- ranger(interactions ~ ., data=web1.df.extended[, -c(1, 2)], importance="impurity")
rf2
sort(importance(rf2), decreasing=T)


###################################################
### code chunk number 26: plot ranger predictions web2
###################################################
preds2.ranger <- predict(rf2, data=web2.df.extended)$predictions
cor(preds2.ranger, web2.df$interactions)
par(mar=c(5,5,1,1))
plot(preds2.ranger +1 , web2.df$interactions + 1, log="xy", las=1, xlab="predicted number 
     of interactions + 1", ylab="observed number of interactions + 1")
abline(0,1)


###################################################
### code chunk number 27: CaseStudy.Rnw:359-361
###################################################
dmultinom(as.vector(tapnet_web2$networks[[1]]$web), prob=preds2.ranger/sum(preds2.ranger), 
          size=sum(tapnet_web2$networks[[1]]$web), log=T)


###################################################
### code chunk number 28: set up tapnets
###################################################
tapnet_web1 <- make_tapnet(tree_low=plant_tree, tree_high=humm_tree, networks=networks[1],
                           traits_low=plant_traits, traits_high=humm_traits, 
                           abun_low=plant_abun[1], abun_high=humm_abun[1], npems_lat=4)
tapnet_web2 <- make_tapnet(tree_low=plant_tree, tree_high=humm_tree, networks=networks[2],
                           traits_low=plant_traits, traits_high=humm_traits, 
                           abun_low=plant_abun[2], abun_high=humm_abun[2], npems_lat=4)
tapnet_web3 <- make_tapnet(tree_low=plant_tree, tree_high=humm_tree, networks=networks[3],
                           traits_low=plant_traits, traits_high=humm_traits, 
                           abun_low=plant_abun[3], abun_high=humm_abun[3], npems_lat=4)


###################################################
### code chunk number 29: tapnet fits
###################################################
fit_web1 <- fit_tapnet(tapnet = tapnet_web1, method="SANN")
fit_web2 <- fit_tapnet(tapnet = tapnet_web2, method="SANN")
fit_web3 <- fit_tapnet(tapnet = tapnet_web3, method="SANN")


###################################################
### code chunk number 30: tapnet fit ells
###################################################
-c(fit_web1$opt$value, fit_web2$opt$value, fit_web3$opt$value)


###################################################
### code chunk number 31: tapnet predict to others
###################################################
preds2.tapnet1 <- predict_tapnet(fit=fit_web1, abuns=tapnet_web2$networks[[1]]$abuns)
preds3.tapnet1 <- predict_tapnet(fit=fit_web1, abuns=tapnet_web3$networks[[1]]$abuns)
preds1.tapnet2 <- predict_tapnet(fit=fit_web2, abuns=tapnet_web1$networks[[1]]$abuns)
preds3.tapnet2 <- predict_tapnet(fit=fit_web2, abuns=tapnet_web3$networks[[1]]$abuns)
preds1.tapnet3 <- predict_tapnet(fit=fit_web3, abuns=tapnet_web1$networks[[1]]$abuns)
preds2.tapnet3 <- predict_tapnet(fit=fit_web3, abuns=tapnet_web2$networks[[1]]$abuns)


###################################################
### code chunk number 32: tapnet compute cors
###################################################
cors.tapnet <- c(
  cor(as.vector(preds2.tapnet1), as.vector(tapnet_web2$networks[[1]]$web)),
  cor(as.vector(preds3.tapnet1), as.vector(tapnet_web3$networks[[1]]$web)),
  cor(as.vector(preds1.tapnet2), as.vector(tapnet_web1$networks[[1]]$web)),
  cor(as.vector(preds3.tapnet2), as.vector(tapnet_web3$networks[[1]]$web)),
  cor(as.vector(preds1.tapnet3), as.vector(tapnet_web1$networks[[1]]$web)),
  cor(as.vector(preds2.tapnet3), as.vector(tapnet_web2$networks[[1]]$web))
)
cors.tapnet


###################################################
### code chunk number 33: tapnet compute cvlogliks
###################################################
ellCV.tapnet <- c(
  dmultinom(as.vector(tapnet_web2$networks[[1]]$web), prob=as.vector(preds2.tapnet1), 
            size=sum(tapnet_web2$networks[[1]]$web), log=T),
  dmultinom(as.vector(tapnet_web3$networks[[1]]$web), prob=as.vector(preds3.tapnet1), 
            size=sum(tapnet_web3$networks[[1]]$web), log=T),
  dmultinom(as.vector(tapnet_web1$networks[[1]]$web), prob=as.vector(preds1.tapnet2), 
            size=sum(tapnet_web1$networks[[1]]$web), log=T),
  dmultinom(as.vector(tapnet_web3$networks[[1]]$web), prob=as.vector(preds3.tapnet2), 
            size=sum(tapnet_web3$networks[[1]]$web), log=T),
  dmultinom(as.vector(tapnet_web1$networks[[1]]$web), prob=as.vector(preds1.tapnet3), 
            size=sum(tapnet_web1$networks[[1]]$web), log=T),
  dmultinom(as.vector(tapnet_web2$networks[[1]]$web), prob=as.vector(preds2.tapnet3), 
            size=sum(tapnet_web2$networks[[1]]$web), log=T)
)
ellCV.tapnet


###################################################
### code chunk number 34: abunsonly fit
###################################################
preds1.abunonly <- (tapnet_web1$networks[[1]]$abuns$low / 
        sum(tapnet_web1$networks[[1]]$abuns$low)) %*% t(tapnet_web1$networks[[1]]$abuns$high / 
        sum(tapnet_web1$networks[[1]]$abuns$high)) / sum(tapnet_web1$networks[[1]]$web)
preds2.abunonly <- (tapnet_web2$networks[[1]]$abuns$low / 
      sum(tapnet_web2$networks[[1]]$abuns$low)) %*% t(tapnet_web2$networks[[1]]$abuns$high / 
      sum(tapnet_web2$networks[[1]]$abuns$high)) / sum(tapnet_web2$networks[[1]]$web)
preds3.abunonly <- (tapnet_web3$networks[[1]]$abuns$low / 
      sum(tapnet_web3$networks[[1]]$abuns$low)) %*% t(tapnet_web3$networks[[1]]$abuns$high / 
      sum(tapnet_web3$networks[[1]]$abuns$high)) / sum(tapnet_web3$networks[[1]]$web)


###################################################
### code chunk number 35: abunsonly cors
###################################################
cors.abun <- c(
  cor(as.vector(preds2.abunonly), as.vector(tapnet_web2$networks[[1]]$web)),
  cor(as.vector(preds3.abunonly), as.vector(tapnet_web3$networks[[1]]$web)),
  cor(as.vector(preds1.abunonly), as.vector(tapnet_web1$networks[[1]]$web)),
  cor(as.vector(preds3.abunonly), as.vector(tapnet_web3$networks[[1]]$web)),
  cor(as.vector(preds1.abunonly), as.vector(tapnet_web1$networks[[1]]$web)),
  cor(as.vector(preds2.abunonly), as.vector(tapnet_web2$networks[[1]]$web))
)
cors.abun


###################################################
### code chunk number 36: abunsonly ellCV
###################################################
ellCV.abuns <- c(
  dmultinom(as.vector(tapnet_web2$networks[[1]]$web), prob=as.vector(preds2.abunonly), 
            size=sum(tapnet_web2$networks[[1]]$web), log=T),
  dmultinom(as.vector(tapnet_web3$networks[[1]]$web), prob=as.vector(preds3.abunonly), 
            size=sum(tapnet_web3$networks[[1]]$web), log=T),
  dmultinom(as.vector(tapnet_web1$networks[[1]]$web), prob=as.vector(preds1.abunonly), 
            size=sum(tapnet_web1$networks[[1]]$web), log=T),
  dmultinom(as.vector(tapnet_web3$networks[[1]]$web), prob=as.vector(preds3.abunonly), 
            size=sum(tapnet_web3$networks[[1]]$web), log=T),
  dmultinom(as.vector(tapnet_web1$networks[[1]]$web), prob=as.vector(preds1.abunonly), 
            size=sum(tapnet_web1$networks[[1]]$web), log=T),
  dmultinom(as.vector(tapnet_web2$networks[[1]]$web), prob=as.vector(preds2.abunonly), 
            size=sum(tapnet_web2$networks[[1]]$web), log=T)
)
ellCV.abuns


###################################################
### code chunk number 37: CaseStudy.Rnw:473-476
###################################################
web3.df <- tapnet2df(tapnet_web3)
web3.df.extended <- cbind.data.frame(web3.df, "match"=(web3.df$traitHBill_length_mean_mm - 
                                                web3.df$traitLCorolla_length_mm)^2 )


###################################################
### code chunk number 38: CaseStudy.Rnw:479-491
###################################################
gam1 <- gam(interactions ~ s(pemLV_1, pemHV_1, bs="ts", k=24) +
  s(pemLV_2, pemHV_3, bs="ts", k=24) + s(traitLCorolla_length_mm, k=3) + 
  s(traitHBill_length_mean_mm, k=3) + s(match, k=3) + s(abunL, k=3) + 
    s(abunH, k=3), data=web1.df.extended, family=nb, gamma=1.4)
gam2 <- gam(interactions ~ s(pemLV_1, pemHV_1, bs="ts", k=24) +
  s(pemLV_2, pemHV_3, bs="ts", k=24) + s(traitLCorolla_length_mm, k=3) + 
  s(traitHBill_length_mean_mm, k=3) + s(match, k=3) + s(abunL, k=3) + 
  s(abunH, k=3), data=web2.df.extended, family=nb, gamma=1.4)
gam3 <- gam(interactions ~ s(pemLV_1, pemHV_1, bs="ts", k=24) +
  s(pemLV_2, pemHV_3, bs="ts", k=24) + s(traitLCorolla_length_mm, k=3) + 
  s(traitHBill_length_mean_mm, k=3) + s(match, k=3) + s(abunL, k=3) + 
  s(abunH, k=3), data=web3.df.extended, family=nb, gamma=1.4)


###################################################
### code chunk number 39: CaseStudy.Rnw:493-499
###################################################
preds2.gam1 <- predict(gam1, newdata=web2.df.extended, type="response")
preds3.gam1 <- predict(gam1, newdata=web3.df.extended, type="response")
preds1.gam2 <- predict(gam2, newdata=web1.df.extended, type="response")
preds3.gam2 <- predict(gam2, newdata=web3.df.extended, type="response")
preds1.gam3 <- predict(gam3, newdata=web1.df.extended, type="response")
preds2.gam3 <- predict(gam3, newdata=web2.df.extended, type="response")


###################################################
### code chunk number 40: CaseStudy.Rnw:501-509
###################################################
cors.gam <- c(
  cor(preds2.gam1, web2.df$interactions),
  cor(preds3.gam1, web3.df$interactions),
  cor(preds1.gam2, web1.df$interactions),
  cor(preds3.gam2, web3.df$interactions),
  cor(preds1.gam3, web1.df$interactions),
  cor(preds2.gam3, web2.df$interactions)
)


###################################################
### code chunk number 41: CaseStudy.Rnw:511-525
###################################################
ellCV.gam <- c(
  dmultinom(as.vector(tapnet_web2$networks[[1]]$web), prob=preds2.gam1/sum(preds2.gam1), 
            size=sum(tapnet_web2$networks[[1]]$web), log=T),
  dmultinom(as.vector(tapnet_web3$networks[[1]]$web), prob=preds3.gam1/sum(preds3.gam1), 
            size=sum(tapnet_web3$networks[[1]]$web), log=T),
  dmultinom(as.vector(tapnet_web1$networks[[1]]$web), prob=preds1.gam2/sum(preds1.gam2), 
            size=sum(tapnet_web1$networks[[1]]$web), log=T),
  dmultinom(as.vector(tapnet_web3$networks[[1]]$web), prob=preds3.gam2/sum(preds3.gam2), 
            size=sum(tapnet_web3$networks[[1]]$web), log=T),
  dmultinom(as.vector(tapnet_web1$networks[[1]]$web), prob=preds1.gam3/sum(preds1.gam3), 
            size=sum(tapnet_web1$networks[[1]]$web), log=T),
  dmultinom(as.vector(tapnet_web2$networks[[1]]$web), prob=preds2.gam3/sum(preds2.gam3), 
            size=sum(tapnet_web2$networks[[1]]$web), log=T)
)


###################################################
### code chunk number 42: CaseStudy.Rnw:527-529
###################################################
cors.gam
ellCV.gam


###################################################
### code chunk number 43: CaseStudy.Rnw:536-545
###################################################
tapnet_web1 <- make_tapnet(tree_low=plant_tree, tree_high=humm_tree, networks=networks[1],
              traits_low=plant_traits, traits_high=humm_traits, abun_low=plant_abun[1], 
              abun_high=humm_abun[1], npems_lat=NULL, use.all.pems=T)
tapnet_web2 <- make_tapnet(tree_low=plant_tree, tree_high=humm_tree, networks=networks[2],
              traits_low=plant_traits, traits_high=humm_traits, abun_low=plant_abun[2], 
              abun_high=humm_abun[2], npems_lat=NULL, use.all.pems=T)
tapnet_web3 <- make_tapnet(tree_low=plant_tree, tree_high=humm_tree, networks=networks[3],
              traits_low=plant_traits, traits_high=humm_traits, abun_low=plant_abun[3], 
              abun_high=humm_abun[3], npems_lat=NULL, use.all.pems=T)


###################################################
### code chunk number 44: CaseStudy.Rnw:547-556
###################################################
web1.df <- tapnet2df(tapnet_web1)
web1.df.extended <- cbind.data.frame(web1.df, "match"=(web1.df$traitHBill_length_mean_mm - 
                                                  web1.df$traitLCorolla_length_mm)^2 )
web2.df <- tapnet2df(tapnet_web2)
web2.df.extended <- cbind.data.frame(web2.df, "match"=(web2.df$traitHBill_length_mean_mm - 
                                                  web2.df$traitLCorolla_length_mm)^2 )
web3.df <- tapnet2df(tapnet_web3)
web3.df.extended <- cbind.data.frame(web3.df, "match"=(web3.df$traitHBill_length_mean_mm - 
                                                  web3.df$traitLCorolla_length_mm)^2 )


###################################################
### code chunk number 45: CaseStudy.Rnw:559-562
###################################################
rf1 <- ranger(interactions ~ ., data=web1.df.extended[, -c(1, 2)], importance="impurity")
rf2 <- ranger(interactions ~ ., data=web2.df.extended[, -c(1, 2)], importance="impurity")
rf3 <- ranger(interactions ~ ., data=web3.df.extended[, -c(1, 2)], importance="impurity")


###################################################
### code chunk number 46: CaseStudy.Rnw:564-567
###################################################
head(sort(round(importance(rf1)), decreasing=T))
head(sort(round(importance(rf2)), decreasing=T))
head(sort(round(importance(rf3)), decreasing=T))


###################################################
### code chunk number 47: CaseStudy.Rnw:569-575
###################################################
preds2.ranger1 <- predict(rf1, data=web2.df.extended)$predictions
preds3.ranger1 <- predict(rf1, data=web3.df.extended)$predictions
preds1.ranger2 <- predict(rf2, data=web1.df.extended)$predictions
preds3.ranger2 <- predict(rf2, data=web3.df.extended)$predictions
preds1.ranger3 <- predict(rf3, data=web1.df.extended)$predictions
preds2.ranger3 <- predict(rf3, data=web2.df.extended)$predictions


###################################################
### code chunk number 48: CaseStudy.Rnw:577-585
###################################################
cors.rf <- c(
  cor(preds2.ranger1, web2.df$interactions),
  cor(preds3.ranger1, web3.df$interactions),
  cor(preds1.ranger2, web1.df$interactions),
  cor(preds3.ranger2, web3.df$interactions),
  cor(preds1.ranger3, web1.df$interactions),
  cor(preds2.ranger3, web2.df$interactions)
)


###################################################
### code chunk number 49: CaseStudy.Rnw:587-601
###################################################
ellCV.rf <- c(
  dmultinom(as.vector(tapnet_web2$networks[[1]]$web), prob=preds2.ranger1/sum(preds2.ranger1), 
            size=sum(tapnet_web2$networks[[1]]$web), log=T),
  dmultinom(as.vector(tapnet_web3$networks[[1]]$web), prob=preds3.ranger1/sum(preds3.ranger1), 
            size=sum(tapnet_web3$networks[[1]]$web), log=T),
  dmultinom(as.vector(tapnet_web1$networks[[1]]$web), prob=preds1.ranger2/sum(preds1.ranger2), 
            size=sum(tapnet_web1$networks[[1]]$web), log=T),
  dmultinom(as.vector(tapnet_web3$networks[[1]]$web), prob=preds3.ranger2/sum(preds3.ranger2), 
            size=sum(tapnet_web3$networks[[1]]$web), log=T),
  dmultinom(as.vector(tapnet_web1$networks[[1]]$web), prob=preds1.ranger3/sum(preds1.ranger3), 
            size=sum(tapnet_web1$networks[[1]]$web), log=T),
  dmultinom(as.vector(tapnet_web2$networks[[1]]$web), prob=preds2.ranger3/sum(preds2.ranger3), 
            size=sum(tapnet_web2$networks[[1]]$web), log=T)
)


###################################################
### code chunk number 50: CaseStudy.Rnw:603-605
###################################################
cors.rf
ellCV.rf


###################################################
### code chunk number 51: CaseStudy.Rnw:611-616
###################################################
all.cors.res <- rbind(cors.tapnet, cors.abun, cors.rf, cors.gam)
all.cors.res <- cbind(all.cors.res, rowMeans(all.cors.res))
colnames(all.cors.res) <- c("1 to 2", "1 to 3", "2 to 1", "2 to 3", "3 to 1", "3 to 2", 
                            "average")
round(all.cors.res, 2)


###################################################
### code chunk number 52: CaseStudy.Rnw:619-623
###################################################
all.ellCV.res <- rbind(ellCV.tapnet, ellCV.abuns, ellCV.rf, ellCV.gam)
all.ellCV.res <- cbind(all.ellCV.res, rowMeans(all.ellCV.res))
colnames(all.ellCV.res) <- colnames(all.cors.res)
round(all.ellCV.res)


###################################################
### code chunk number 53: CaseStudy.Rnw:631-659
###################################################
fits.cors.res <- rbind(
  "tapnet"=cbind(
    cor(as.vector(tapnet_web1$networks[[1]]$web), as.vector(predict_tapnet(fit_web1, 
                                           abuns=tapnet_web1$networks[[1]]$abuns))),
    cor(as.vector(tapnet_web2$networks[[1]]$web), as.vector(predict_tapnet(fit_web2, 
                                           abuns=tapnet_web2$networks[[1]]$abuns))),
    cor(as.vector(tapnet_web3$networks[[1]]$web), as.vector(predict_tapnet(fit_web3, 
                                           abuns=tapnet_web3$networks[[1]]$abuns)))
  ),
  "abuns"=cbind(
        cor(as.vector(tapnet_web1$networks[[1]]$web), as.vector(preds1.abunonly)),
        cor(as.vector(tapnet_web2$networks[[1]]$web), as.vector(preds2.abunonly)),
        cor(as.vector(tapnet_web3$networks[[1]]$web), as.vector(preds3.abunonly))
  ),
  "rf"=cbind(
    cor(predict(rf1, data=web1.df.extended)$predictions, web1.df$interactions),
    cor(predict(rf2, data=web2.df.extended)$predictions, web2.df$interactions),
    cor(predict(rf3, data=web3.df.extended)$predictions, web3.df$interactions)
  ),
  "gam"=cbind(
    cor(predict(gam1, data=web1.df.extended, type="response"), web1.df$interactions),
    cor(predict(gam2, data=web2.df.extended, type="response"), web2.df$interactions),
    cor(predict(gam3, data=web3.df.extended, type="response"), web3.df$interactions)
  )
)
fits.cors.res <- cbind(fits.cors.res, rowMeans(fits.cors.res))
colnames(fits.cors.res) <- c("1 to 1", "2 to 2", "3 to 3", "average")
round(fits.cors.res, 2)


###################################################
### code chunk number 54: tapnet without abundances
###################################################
tapnet_web1.w <- make_tapnet(tree_low=plant_tree, tree_high=humm_tree, networks=
              networks[1], traits_low=plant_traits, traits_high=humm_traits, npems_lat=4)
tapnet_web2.w <- make_tapnet(tree_low=plant_tree, tree_high=humm_tree, networks=
              networks[2], traits_low=plant_traits, traits_high=humm_traits, npems_lat=4)
tapnet_web3.w <- make_tapnet(tree_low=plant_tree, tree_high=humm_tree, networks=
              networks[3], traits_low=plant_traits, traits_high=humm_traits, npems_lat=4)


###################################################
### code chunk number 55: fit without abundances
###################################################
fit_web1.w <- fit_tapnet(tapnet=tapnet_web1.w)
fit_web2.w <- fit_tapnet(tapnet=tapnet_web2.w)
fit_web3.w <- fit_tapnet(tapnet=tapnet_web3.w)


###################################################
### code chunk number 56: predict without abundances
###################################################
preds2.tapnet1.w <- predict_tapnet(fit=fit_web1.w, abuns=tapnet_web2.w$networks[[1]]$abuns)
preds3.tapnet1.w <- predict_tapnet(fit=fit_web1.w, abuns=tapnet_web3.w$networks[[1]]$abuns)
preds1.tapnet2.w <- predict_tapnet(fit=fit_web2.w, abuns=tapnet_web1.w$networks[[1]]$abuns)
preds3.tapnet2.w <- predict_tapnet(fit=fit_web2.w, abuns=tapnet_web3.w$networks[[1]]$abuns)
preds1.tapnet3.w <- predict_tapnet(fit=fit_web3.w, abuns=tapnet_web1.w$networks[[1]]$abuns)
preds2.tapnet3.w <- predict_tapnet(fit=fit_web3.w, abuns=tapnet_web2.w$networks[[1]]$abuns)


###################################################
### code chunk number 57: tapnet without abundances compute cors
###################################################
cors.tapnet.w <- c(
  cor(as.vector(preds2.tapnet1.w), as.vector(tapnet_web2.w$networks[[1]]$web)),
  cor(as.vector(preds3.tapnet1.w), as.vector(tapnet_web3.w$networks[[1]]$web)),
  cor(as.vector(preds1.tapnet2.w), as.vector(tapnet_web1.w$networks[[1]]$web)),
  cor(as.vector(preds3.tapnet2.w), as.vector(tapnet_web3.w$networks[[1]]$web)),
  cor(as.vector(preds1.tapnet3.w), as.vector(tapnet_web1.w$networks[[1]]$web)),
  cor(as.vector(preds2.tapnet3.w), as.vector(tapnet_web2.w$networks[[1]]$web))
)
cors.tapnet.w
mean(cors.tapnet.w)


