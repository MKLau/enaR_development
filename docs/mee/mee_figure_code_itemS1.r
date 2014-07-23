# enaR Example Code
# Generates enaR example plot Figure for MEE paper
# Borrett & Lau
# July 2013
# -----------------------------------------------

# preliminaires
rm(list=ls())  # clear working memory
library(enaR)  # load enaR package  (this script usese version 2.0)
data(troModels)  # load trophic model library
data(oyster)     # load oyster model separately for convenince

# start figure
#fn <- "../figures/enaR_plot_example.eps"
fn <- "enaR_plot_example.eps"
postscript(fn,width=6,height=5.5,
           onefile=TRUE,
           horizontal=FALSE,
           paper="special",
           family="Times",pointsize=12)
#pdf(fn,width=6,height=5.5,family="Times",pointsize=12)


opar <- par(las=1,mfrow=c(2,2),mar=c(4,5,1,1),oma=c(1,1,1,1))
# ----------------------------------
# panel a -- network visualization
# ----------------------------------
my.col=c("red4","grey", "black")
F=oyster%n%'flow'                   # extract flow information for later use.
f=which(F!=0, arr.ind=T)       # get indices of positive flows
opar <- par(xpd=TRUE)  #,mai=c(1.02, 0.62, 0.82, 0.42))
set.seed(2)                    # each time the plot is called, the
                               # layout orientation changes.  setting
                               # the seed ensures a consistent
                               # orientation each time the plot
                               # function is called.
# 
plot(oyster,
     vertex.cex=log(oyster%v%'storage'), # scale nodes with storage
     label= oyster%v%'vertex.names',     # add node labels
     boxed.labels=FALSE,
     label.cex=0.8,     # scale label text
     vertex.sides=45,   # to make rounded
     edge.lwd=log10(abs(F[f]))+2,     # scale arrows to flow magnitude
     edge.col=my.col[3],  # change the edge color
     vertex.col=my.col[1],   # change the vertx color
     label.col="black",      # set teh label color
     vertex.border = my.col[3],
     vertex.lty = 1,         # vertex line type
     xlim=c(-4,1),ylim=c(-2,-2))
#
mtext("(a)",side=3,line=0,adj=-0.5)
opar <- par(xpd=FALSE)
#
# ------------------------
# panel b -- HMG bargraph
# ------------------------
#
# balance models
m.list <- lapply(troModels,force.balance)
#
# batch apply flow analysis to 56 models and extract HMG network statistic
hmg <- unlist(lapply(m.list,function(x) enaFlow(x)$ns[13]))
hmg <- sort(hmg,decreasing=TRUE) # sort hmg values
#
# create bar plot
barplot(hmg,
        ylim=c(0,2),
        col="darkgreen",border=NA,
#        xlab="Trophic Models (Rank Ordered)",
#        ylab="Network Homogenization",
        ylab="",
        names.arg="")
mtext("Trophic Ecosystem Models (ranked)",side=1,line=1,cex=1)
mtext("Network Homogenization",side=2,line=2.5,cex=1,las=3)
abline(h=1,col="orange",lwd=2)
mtext("(b)",side=3,line=0,adj=-0.5)
#
# ------------------------
# panel c -- I/D vs A/C
# ------------------------
#
# Collecting and combining all network statistics
ns.list <- lapply(m.list,get.ns) # returns as list
ns <- do.call(rbind,ns.list)  # ns as a data.frame
#
plot(ns$ASC.CAP,ns$IFI,pch=20,col="blue",cex=2,
     ylab="",
     xlab="",
     xlim=c(0,0.8),ylim=c(0,0.8))
mtext("Indirect Flow Intensity (IFI)",side=2,line=2.5,cex=1,las=3)
mtext("Ascendency-to-Capacity (A/C)",side=1,line=2.5,cex=1)
mtext("(c)",side=3,line=0,adj=-0.5)
#
# ----------------------------------
# panel d -- Centrality Target Plot
# ----------------------------------
#
m <- troModels[[38]]  # extract Chesapeake Bay Model
b <- betweenness(m)         # calculate betweenness centrality (SNA)
nms <- m%v%'vertex.names'   # get vertex names
nms[35] <- "Susp. Part. Carbon"
nms[36] <- "Sed. Part. Carbon"
nms[b<=(0.1*max(b))] <- NA  # exclude less central nodes
#
set.seed(3)
opar <- par(xpd=TRUE,mar=c(3,5,1,3))
# create target plot
gplot.target(m,b,#circ.lab=FALSE,
             edge.col="grey",
             vertex.cex=1.1,
             label=nms, # show only labels of most central nodes
             circ.lab=FALSE,
             label.cex=0.75)
mtext("(d)",side=3,line=0,adj=-0.5)
                                        #xlim=c(-1,4))
# ---
dev.off()  # close eps file
cmd <- paste("open ",fn)
system(cmd) # open eps file (on mac)

