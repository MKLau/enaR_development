###enaR function output checking algorithm
###Compares the output of for the Cone Springs Model
###MKLau 23Jul2014

library(devtools)
install_github(username='TheSeeLab',repo='SEELab/enaR',ref='beta')
library(enaR)
data(enaModels)
enaModels <- lapply(enaModels,force.balance)
ns.beta <- lapply(enaModels,get.ns)
dput(do.call(rbind,ns.beta),file='~/projects/packages/enaR_development/data/nsbeta.rda')

q('no')
library(devtools)
install_github(username='TheSeeLab',repo='SEELab/enaR',ref='edgeflow')
library(enaR)
data(enaModels)
enaModels <- lapply(enaModels,force.balance)
ns.edgeflow <- lapply(enaModels,get.ns)
dput(do.call(rbind,ns.edgeflow),file='~/projects/packages/enaR_development/data/nsedgeflow.rda')

###Checks
dget(file='~/projects/packages/enaR_development/data/nsbeta.rda')
dget(file='~/projects/packages/enaR_development/data/nsedgeflow.rda')
