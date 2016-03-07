###Adding more model info to enaModelInfo
###mklau 29 Sep 2014

library(enaR)
data(enaModels)

###Loading
emi <- read.csv('~/projects/packages/enaR_development/data/enaModelInfo.csv')
emi <- emi[emi$enaR==1,]
emi[,1] <- str_replace(pattern='\'\'',replacement='',string=emi[,1])
emi[,1] <- str_replace(pattern='\'\'',replacement='',string=emi[,1])
emi[,1] <- as.character(emi[,1])
emi <- emi[is.na(emi$Model)==FALSE,]

###Coallesce with the model set
emi$Model <- sort(names(enaModels))
all(emi$Model==names(enaModels))

###Sort
emi <- emi[match(names(enaModels),emi$Model),]
all(emi$Model==names(enaModels))

###De-limit
emi <- emi[,1:8]

###Save
detach(package:enaR)
enaModelInfo <- emi
package.skeleton(path='~/projects/packages/enaR_development/data/',name='emi')
