library(enaR)
library(intergraph)
library(igraph)
library(Rgraphviz)

enaPlot <- function(x,v.col='white',e.col='grey'){
    g1 <- asIgraph(x)
    x <- unpack(x)
    V(g1)$name <- rownames(x$F)
    g1 <- as_graphnel(g1)
    glty <- x$living
    if (all(is.na(glty))){glty <- rep(1,length(glty))}else{
        glty <- sign(glty)
        glty[glty == 0] <- 2
    }
    names(glty) <- rownames(x$F)
    gl1 <- layoutGraph(g1)
    nodeRenderInfo(gl1) <- list(lty=glty)
    attrs <- list(node=list(fillcolor=v.col),edge=list(color=e.col))
    plot(gl1,attrs=attrs)
}

