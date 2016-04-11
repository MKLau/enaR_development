### Entering the pitcher plant model from Mouquet 2008
### MKLau 11Apr2016

### Functional Ecology, August, 2008. 10.1111/j.1365-2435.2008.01421.x

library(enaR)
library(intergraph)
library(igraph)
library(Rgraphviz)
source('../src/helpers.R')

### model variables
                                        #inputs (mg C L^-1 Day^-1)
i <- list(A=5.39)
                                        #consumption values
u <- list(B=0.001,P=0.014,M=0.5872)
                                        #mortality values (day^-1)
m <- list(B=0.001,P=0.01)
                                        #respiration values
r <- list(B=0.0005,P=0.0014,M=0.01)
                                        #sedimentation values (day^-1)
s <- list(D=0.01)
                                        #export (day^-1)
e <- list(M=1)

### constructing the flow matrix
                                        # edgelist
E <- data.frame(matrix(c(
    'ants','detritus',1,
    'detritus','sediment',s$D,
    'detritus','bacteria',u$B,
    'bacteria','protozoa',u$P/2,
    'bacteria','rotifers',u$P/2,
    'rotifers','mosquitos',u$M/2,
    'protozoa','mosquitos',u$M/2,
    'bacteria','detritus',m$B,
    'protozoa','detritus',m$P/2,
    'rotifers','detritus',m$P/2,
    'bacteria','respiration',r$B,
    'protozoa','respiration',r$P/2,
    'rotifers','respiration',r$P/2,
    'mosquitos','respiration',r$M,
    'mosquitos','export',e$M
    ),ncol=3,byrow=TRUE)
)
E[,3] <- as.numeric(as.character(E[,3]))

flow <- edge2flow(E)

resp <- flow[-(nrow(flow)-1):-nrow(flow),colnames(flow)=='respiration']
expo <- flow[-(nrow(flow)-1):-nrow(flow),colnames(flow)=='export']
flow.m <- flow[-(nrow(flow)-1):-nrow(flow),
                 -(ncol(flow)-1):-ncol(flow)]
inpu <- c(i$A,rep(0,nrow(flow.m)-1))
stor <- c(A=0,D=371.25,B=5.00,R=27.85/2,P=27.85/2,M=303.75,S=0)
livi <- c(A=TRUE,D=FALSE,B=TRUE,R=TRUE,P=TRUE,M=TRUE,S=FALSE)
moq.mod <- pack(f=flow.m,inp=inpu,res=resp,exp=expo,sto=stor,liv=livi)

### 

