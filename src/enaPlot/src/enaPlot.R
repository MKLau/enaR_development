library(enaR)
library(intergraph)
library(igraph)
library(Rgraphviz)

enaPlot <- function(x,v.col='white',e.col='grey',extended=FALSE,lab.size.limit=FALSE,max.lab=5){
    if (!(lab.size.limit)){max.lab <- max(nchar(rownames(unpack(x)$F)))}
    if (extended){
        x. <- unpack(x)
        I <- x.$z
        E <- R <- rep(0,nrow(x.$F))
        ef <- rbind(x.$F,E,R,I)
        I <- rep(0,nrow(ef))
        E <- c(x.$e,rep(0,3));R <- c(x.$r,rep(0,3))
        ef <- cbind(ef,E,R,I)
        x <- pack(ef,inp=x.$z,res=x.$r,exp=x.$e,sto=x.$X,liv=x.$living)
    }else{}
    if (lab.size.limit){
        v.labs <- rownames(unpack(x)$F)
        for (i in 1:length(v.labs)){
            if (max.lab < nchar(v.labs[i])){v.labs[i] <- substr(v.labs[i],1,max.lab)}
        }
        if (any(duplicated(v.labs))){
            for (i in 1:length(table(v.labs)[table(v.labs) > 1])){
                v.labs[v.labs == names(table(v.labs))[i]] <- 
                    paste(names(table(v.labs))[i],1:(table(v.labs)[i]),sep='')
            }
        }
        x%v%'vertex.names' <- v.labs
    }else{}
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

