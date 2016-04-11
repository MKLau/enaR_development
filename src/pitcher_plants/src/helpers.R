edge2flow <- function(x){
    x[,1:2] <- apply(x[,1:2],2,as.character)
    x[,3] <- as.numeric(as.character(x[,3]))
    e.l <- unique(unlist(x[,1:2]))
    out <- matrix(0,nrow=length(e.l),ncol=length(e.l))
    rownames(out) <- colnames(out) <- e.l
    for (i in 1:nrow(x)){
        out[rownames(out) == x[i,1],colnames(out) == x[i,2]] <- x[i,3]
    }
    return(out)
}
