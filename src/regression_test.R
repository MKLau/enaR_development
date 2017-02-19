########################################################################################################################
### Remove packages and check which functions break
########################################################################################################################


## Depends:
##     R (>= 2.10),
##     MASS,
##     stringr,
##     sna,
##     network,
##     gdata,
##     stats,
##     utils,
##     graphics,
##     limSolve
## Suggests:
##     codetools,
##     igraph,
##     R.rsp
## VignetteBuilder: R.rsp

library(enaR)
data(enaModels)

ls('package:enaR', all.names=TRUE)
sessionInfo()

ShannonDiversity(runif(10))
TES(enaModels[[1]])
TET(enaModels[[1]])
as.bipartite()
as.extended()
enaR:::bal()
balance(enaModels[[1]])
enaR:::bcratio(as.matrix(enaModels[[1]]))
enaR:::cycliv(enaModels[[1]])
eigenCentrality(as.matrix(enaModels[[1]]))
enaAll(enaModels[[1]])
enaAscendency(enaModels[[1]])
enaControl(enaModels[[1]])
enaCycle(enaModels[[1]])
enaEnviron(enaModels[[1]])
enaFlow(enaModels[[1]])
enaMTI(enaModels[[1]])
enaStorage(enaModels[[1]])
enaStructure(enaModels[[1]])
enaTroAgg(enaModels[[1]])
enaUtility(enaModels[[1]])
environCentrality(enaModels[[1]])
findPathLength(enaModels[[1]])
balance(enaModels[[1]])
get.ns(enaModels[[1]])
get.orient()
mExp(as.matrix(enaModels[[1]]))
netOrder(enaModels[[1]])
scc(as.matrix(enaModels[[1]]))
get.orient()
enaR:::signs(as.matrix(enaModels[[1]]))
enaR:::structure.statistics(as.matrix(enaModels[[1]]))
ssCheck(enaModels[[1]])
unpack(enaModels[[1]])

## read.enam()
## read.nea()
## read.scor()
## read.wand()
## enaR:::relationalChange()
## enaR:::scifix()
## write.nea()

