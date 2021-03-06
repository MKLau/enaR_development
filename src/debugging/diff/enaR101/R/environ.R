# environ --- conducts environ analysis 
# INPUT = network object
# OUTPUT = input and/or output environs
# 
# M. Lau | July 2011
# ---------------------------------------------------

environ <- function(x = 'network object', input = TRUE, output = TRUE, err.tol = 1e-10,balance.override=FALSE){
                                        #check for network class
  if (class(x) != 'network'){warning('x is not a network class object')}

                                          #Check for balancing
  if (balance.override){}else{
    if (any(list.network.attributes(x) == 'balanced') == FALSE){x%n%'balanced' <- ssCheck(x)}
    if (x%n%'balanced' == FALSE){warning('Model is not balanced'); stop}
  }

  F <- enaFlow(x)   # calculate enaFlow

                                        #calculations
  if (input == TRUE){
                                        #Input perspective
    EP <- list()
    for (i in (1:nrow(F$NP))){
      dNP <- diag(F$NP[i,]) #diagonalized N matrix
      EP[[i]] <- dNP %*% F$GP         # calculate internal environ flows
      EP[[i]] <- EP[[i]] - dNP      # place negative environ throughflows on the principle diagonal
      EP[[i]] <- cbind(EP[[i]],apply((-EP[[i]]),1,sum))
      EP[[i]] <- rbind(EP[[i]],c(apply((-EP[[i]]),2,sum)))
      EP[[i]][nrow(EP[[i]]),ncol(EP[[i]])] <- 0
      EP[[i]][abs(EP[[i]]) < err.tol] <- 0
    }
    names(EP) <- rownames(F$GP)
  }
  if (output == TRUE){
                                        #Output perspective
    E <- list()
    for (i in (1:nrow(F$N))){
      dN <- diag(F$N[,i]) #diagonalized N matrix
      E[[i]] <- F$G %*% dN
      E[[i]] <- E[[i]] - dN
      E[[i]] <- cbind(E[[i]],apply((-E[[i]]),1,sum))
      E[[i]] <- rbind(E[[i]],c(apply((-E[[i]]),2,sum)))
      E[[i]][nrow(E[[i]]),ncol(E[[i]])] <- 0
      E[[i]][abs(E[[i]]) < err.tol] <- 0
    }
    names(E) <- rownames(F$G)
  }
                                        #Wrap-up output into list  
  if (input == TRUE & output == TRUE){
    out <- list('input' = EP,'output' = E)
  }else if (input == TRUE & output == FALSE){
    out <- EP
  }else if (input == FALSE & output == TRUE){
    out <- E
  }

  return(out)

}
