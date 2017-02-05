### Reading models from the econet website
### MKLau 31 Mar 2016







#' Access example EcoNet models from the EcoNet website.
#' 
#' This function allows the user to access models that are presented on the
#' website for EcoNet, the web-based interface for conducting ENA
#' (http://eco.engr.uga.edu/), by Caner Kazanci at the University of Georgia.
#' 
#' 
#' @param x The URL for the EcoNet examples.
#' @param model.name The model to be accessed. If 'prompt' the user will be
#' asked for the model they wish to use. Can also be a number for the model or
#' the name of the model.
#' @return Returns the model formatted as a network object.
#' @author Matthew K. Lau
#' @seealso \code{\link{read.EcoNet}}
#' @references Kazanci, C., 2007. EcoNet: A new software for ecological
#' modeling, simulation and network analysis, Ecol. Model., Vol 208/1 pp 3-8.
#' @export EcoNetWeb
EcoNetWeb <- function(x='http://eco.engr.uga.edu/Examples/examples.html',model.name='prompt'){
    x <- readLines(x)
    x <- x[1:(grep('<!-- Copy paste for new model',x)-1)]
    mod.names <- x[grep('<h2>',x)]
    mod.names <- sub('<h2>','',mod.names);mod.names <- sub('</h2>','',mod.names)
    mod.names <- strsplit(mod.names,split=' ')
    mod.names <- lapply(mod.names,function(x) x[x != ''])
    mod.names <- unlist(lapply(mod.names,paste,collapse=' '))
    mod.loc <- cbind(grep('<hr class=\"duz\">',x),grep('</pre></div></td></tr></table>',x))
    mods <- apply(mod.loc,1,function(loc,x) x[loc[1]:loc[2]],x=x)
    if (model.name == 'prompt'){
        print(mod.names)
        model.name <- as.numeric(readline('Enter model number:'))
        out <- read.EcoNet(mods[[model.name]])
    }else if (is.numeric(model.name[1])){
        out <- read.EcoNet(mods[[model.name[1]]])
    }else{
        out <- read.EcoNet(mods[[agrep(model.name[1],mod.names,ignore.case=T)]])
    }
    return(out)
}
#' Shannon Diversity Metrics
#' These are based on entropy and build Shannon and Weaver 1949
#'
#' Borrett | November 29, 2016
#'
#' INPUT = Vector
#' Output = set of network statistics to charcterize the diversity in the vector
#' ================================================================================







#' Shannon Diversity Metrics These are based on entropy and build Shannon and
#' Weaver 1949
#' 
#' Borrett | November 29, 2016
#' 
#' INPUT = Vector Output = set of network statistics to charcterize the
#' diversity in the vector
#' ================================================================================
#' Shannon Diversity Metrics These are based on entropy and build Shannon and
#' Weaver 1949
#' 
#' Borrett | November 29, 2016
#' 
#' INPUT = Vector Output = set of network statistics to charcterize the
#' diversity in the vector
#' ================================================================================
#' Shannon information entropy
#' 
#' Calculates a number of metrics based on the Shannon information entropy
#' measure of diversity in a vector, x.
#' 
#' @param x a 1 x n vector.
#' @return \item{H}{Shannon entropy-based metric of diversity.  This captures
#' the effects of both richnes (the length of the vector, n) and the evenennes
#' of the distribution.} \item{Hmax}{The maximum possible value of H given a
#' vector of the length n provided.} \item{Hr}{Relative evenness Hr = H/Hmax}
#' \item{Hcentral}{The centralization or concentration of the values among the
#' n elements} \item{n}{Number of elements in the vector.}
#' \item{effective.n}{effective number of elements in the vector, given the
#' distribution of the relative weights.}
#' @note The formulation for Shannon Diversity uses a natural logarithm.  As
#' the natural logarithm of zero is undefined, the input vector cannot contain
#' zeros.  Analytically, there are two approaches to dealing with this issue if
#' your vector contains zeros.  First, you can apply the analysis to only the
#' non-zero elements.  Second, you can add a tiny amount to all of the elements
#' such that the zero elements are now very small numbers, relative the
#' original vector values.
#' @author Stuart R. Borrett
#' @export ShannonDiversity
#' @examples
#' 
#' 
#' 
#' data(oyster)
#' 
#' #' throughflow diversity
#' T <- enaFlow(oyster)$T
#' ShannonDiversity(T)
#' 
#' #' storage (biomass) biodiversity
#' X <- oyster %v% "storage"
#' ShannonDiversity(X)
#' 
#' 
#' 
ShannonDiversity <-  function(x){
    p <- x/sum(x)  # relative proportion
    H <- -1 * sum(p * log(p) )  # results in nats (using natural log)  # Shannon Diversity
    Hmax <- log(length(x)) # maximum possible Shannon Diversity
    Hr <- H/Hmax
    Hcentral <- 1-Hr
    effective.n <- exp(H)  # effecive number of elements
    n <- length(x)  # number of elements

    return(c("H"=H, "Hmax" = Hmax, "Hr" = Hr, "Hcentral" = Hcentral,
            "n" = n,
            "effective.n" = effective.n))
}
#' TES.R  --- TOTAL ENVIRON STORAGE
#' INPUT = network model
#' OUTPUT = total environ throughput - unit and scaled
#'
#' Borrett | July 7, 2012
#' ---------------------------------------------------







#' TES.R --- TOTAL ENVIRON STORAGE INPUT = network model OUTPUT = total environ
#' throughput - unit and scaled
#' 
#' Borrett | July 7, 2012 ---------------------------------------------------
#' TES.R --- TOTAL ENVIRON STORAGE INPUT = network model OUTPUT = total environ
#' throughput - unit and scaled
#' 
#' Borrett | July 7, 2012 ---------------------------------------------------
#' Calculate the Total Environ Storage
#' 
#' Calculates the total storage in each n input and output environs.  This
#' function calculates the storage for both the unit input (output) and the
#' realized input (output) environs.  Realized uses the observed inputs
#' (outputs) rather than an assumed unit input (output) to each node.
#' 
#' @param x A network object.
#' @param balance.override LOGICAL: should balancing being ignored.
#' @return \item{realized.input}{input oriented, realized storage in each
#' environ.} \item{realized.output}{output oriented, realized storage in each
#' environ.} \item{unit.input }{input oriented, unit storage in each environ.}
#' \item{unit.output}{input oriented, unit storage in each environ.}
#' @author Matthew K. Lau Stuart R. Borrett David E. Hines
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' \code{\link{enaStorage},\link{enaEnviron}}
#' @references Matis, J.H. and Patten, B.C. 1981.  Environ analysis of linear
#' compartmenal systems: the static, time invariant case.  Bulletin of the
#' International Statistical Institute. 48, 527--565.
#' @examples
#' 
#' 
#' 
#' data(troModels)
#' tes <- TES(troModels[[6]])
#' tes
#' 
#' 
#' 
#' @export TES
TES <- function(x,balance.override=FALSE){

                                        #Check for network class
  if (class(x) != 'network'){warning('x is not a network class object')}

                                        #Check for balancing
  if (balance.override){}else{
    if (any(list.network.attributes(x) == 'balanced') == FALSE){x%n%'balanced' = ssCheck(x)}
    if (x%n%'balanced' == FALSE){warning('Model is not balanced'); stop}
  }
                                        #
  oo <- get.orient() #original orientation
  if (oo == 'school'){oo <- 'internal'}
  set.orient('internal')
  S <- enaStorage(x)
  set.orient(oo)
  input <- unpack(x)$z   # get data input elements
  output <- unpack(x)$y  # get data output elements

  # UNIT
  X = S$S 
  unit.output <- apply(X,2,sum)
  X =  S$SP
  unit.input <- apply(X,2,sum)

  # REALIZED
  X = S$S %*% diag(input)
  realized.output <- apply(X,2,sum)
  X =  diag(output) %*% S$SP
  realized.input <- apply(X,2,sum)
  
  return(list("realized.input"=realized.input,"realized.output"=realized.output,"unit.input"=unit.input,"unit.output"=unit.output))
}
#' TET.R  --- TOTAL ENVIRON THROUGHFLOW
#' INPUT = network model
#' OUTPUT = total environ throughput - unit and scaled
#'
#' Borrett | July 7, 2012
#' ---------------------------------------------------







#' TET.R --- TOTAL ENVIRON THROUGHFLOW INPUT = network model OUTPUT = total
#' environ throughput - unit and scaled
#' 
#' Borrett | July 7, 2012 ---------------------------------------------------
#' TET.R --- TOTAL ENVIRON THROUGHFLOW INPUT = network model OUTPUT = total
#' environ throughput - unit and scaled
#' 
#' Borrett | July 7, 2012 ---------------------------------------------------
#' Calculates the Total Environ Throughflow for a Ecosystem Network Model
#' 
#' Determines the total environ throughflow (TET) for each of the 2 x n
#' environs of the selected network model. It returns both the TET calculated
#' from a unit input (output) vector and from the observed or realized input
#' (output) vector.
#' 
#' @param x A network object.
#' @param balance.override Logical: should the function work if the model is
#' not at steady-state?
#' @return \item{realized.input}{vector of the n realized total environ
#' throughflows for the n input oriented environs.}
#' \item{realzied.output}{vector of the n realized total environ throughflows
#' for the n ouptut oriented environs.} \item{unit.input}{vector of the n unit
#' total environ throughflows for the n input oriented environs.}
#' \item{unit.output}{vector of the n unit total environ throughflows for the n
#' output oriented environs.}
#' @author Matthew K. Lau Stuart R. Borrett
#' @seealso \code{\link{enaEnviron}}
#' @references Gattie, D.K., Schramski, J.R., Borrett, S.R., Patten, B.C.,
#' Bata, S.A., and Whipple, S.J. 2006. Indirect effects and distributed control
#' in ecosystems: Network environ analysis of a seven-compartment model of
#' nitrogen flow in the Neuse River Estuary, USA---Steady-state analysis. Ecol.
#' Model. 194:162--177.
#' 
#' Whipple, S.J., Borrett, S.R., Patten, B.C., Gattie, D.K., Schramski, J.R.,
#' and Bata, S.A. 2007.  Indirect effects and distributed control in
#' ecosystems: Comparative network environ analysis of a seven-compartment
#' model of nitrogen flow in the Neuse River Estuary, USA---Time series
#' analysis. Ecol. Model. 206: 1--17.
#' @examples
#' 
#' 
#' 
#' data(troModels)
#' tet <- TET(troModels[[6]])
#' tet
#' 
#' 
#' 
#' @export TET
TET <- function(x,balance.override=FALSE){

                                        #Check for network class
  if (class(x) != 'network'){warning('x is not a network class object')}

                                        #Check for balancing
  if (balance.override){}else{
    if (any(list.network.attributes(x) == 'balanced') == FALSE){x%n%'balanced' = ssCheck(x)}
    if (x%n%'balanced' == FALSE){warning('Model is not balanced'); stop}
  }

  oo <- get.orient() #original orientation
  if (oo == 'school'){oo <- 'internal'}
  set.orient('internal')
  E <- enaEnviron(x)
  set.orient(oo)
  input <- unpack(x)$z   # get data input elements
  output <- unpack(x)$y  # get data output elements

  # UNIT
  unit.input <- 0
  unit.output <- 0

  # REALIZED
  realized.input <- 0 # initialize
  realized.output <- 0 # initialize

  # UNIT & SCALED
  for(i in 1:length(input)){
    realized.input[i] = -sum(diag(E$input[[i]] * output[i]))
    realized.output[i] = -sum(diag(E$output[[i]] * input[i]))
    unit.input[i] = -sum(diag(E$input[[i]]))
    unit.output[i] = -sum(diag(E$output[[i]]))
  }

  return(
         list("realized.input"=realized.input,
              "realized.output"=realized.output,
              "unit.input"=unit.input,
              "unit.output"=unit.output)
         )

}
#' as.bipartite  --- convert a network object to a matrix for 
#' analysis with the bipartite package
#' INPUT = network model
#' OUTPUT = matrix representation
#' M. Lau July 2015
#' ---------------------------------------------------







#' as.bipartite --- convert a network object to a matrix for analysis with the
#' bipartite package INPUT = network model OUTPUT = matrix representation M.
#' Lau July 2015 ---------------------------------------------------
#' as.bipartite --- convert a network object to a matrix for analysis with the
#' bipartite package INPUT = network model OUTPUT = matrix representation M.
#' Lau July 2015 --------------------------------------------------- Create a
#' bipartite network.
#' 
#' Converts a network object (unipartite) to a two-mode (bipartite) network
#' representation.
#' 
#' Bipartite network approaches are often used for analyzing the structure of
#' interactions among species in communities.  Although typically ecosystem
#' networks are handled using a unipartite representation, anlayzing them from
#' a bipartite perspective may be informative. This function provides an easy
#' means for converting to a bipartite representation as long as there is a
#' natural division to categorize species into distinct modes.
#' 
#' @param x A network object.
#' @param y A vector of membership values.
#' @return Returns a matrix with the species of one mode arrayed in rows and
#' the other in columns.
#' @author Matthew K. Lau
#' @examples
#' 
#' 
#' 
#' data(oyster)
#' as.bipartite(oyster, gl(2, 3))
#' 
#' 
#' 
#' @export as.bipartite
as.bipartite <- function(x,y){
    y <- factor(y)
    unpack(x)$F[y == levels(y)[1],y == levels(y)[2]]
}

#' as.extended  --- convert a network object to extended format
#' in Allesina and Bondavalli 2003
#' INPUT = network model
#' OUTPUT = the same model in extended format with inputs and
#' exports/respiration in the same matrix
#' REFERENCE: Allesina, S., Bondavalli, C., 2003.
#' Steady state of ecosystem flow networks: a comparison
#' between balancing procedures. Ecological Modelling 165(2-3):
#' 231-239.
#' M. Lau July 2011
#' ---------------------------------------------------







#' as.extended --- convert a network object to extended format in Allesina and
#' Bondavalli 2003 INPUT = network model OUTPUT = the same model in extended
#' format with inputs and exports/respiration in the same matrix REFERENCE:
#' Allesina, S., Bondavalli, C., 2003. Steady state of ecosystem flow networks:
#' a comparison between balancing procedures. Ecological Modelling 165(2-3):
#' 231-239. M. Lau July 2011
#' --------------------------------------------------- as.extended --- convert
#' a network object to extended format in Allesina and Bondavalli 2003 INPUT =
#' network model OUTPUT = the same model in extended format with inputs and
#' exports/respiration in the same matrix REFERENCE: Allesina, S., Bondavalli,
#' C., 2003. Steady state of ecosystem flow networks: a comparison between
#' balancing procedures. Ecological Modelling 165(2-3): 231-239. M. Lau July
#' 2011 --------------------------------------------------- Create an Extended
#' Format Matrix
#' 
#' Converts a network object to the extended format of Allesina and Bondavalli
#' (2003).
#' 
#' Used in the balance function.
#' 
#' @param x A network object.
#' @param zero.na Logical: should NA's be replaced with zeros?
#' @return Returns an extended format matrix.
#' @author Matthew K. Lau Stuart R. Borrett
#' @seealso \code{\link{balance}}
#' @references Allesina, S., Bondavalli, C., 2003. Steady state of ecosystem
#' flow networks: a comparison between balancing procedures.Ecological
#' Modelling 165(2-3):231-239.
#' @examples
#' 
#' 
#' 
#' data(troModels)
#' as.extended(troModels[[6]])
#' 
#' 
#' 
#' @export as.extended
as.extended <- function(x,zero.na=TRUE){
                                        #Check for network class object
  if (class(x) != "network"){warning('x is not a network class object')}
                                        #unpack the data from the network object
  flow <- as.matrix(x,attrname="flow")
  input <- x%v%'input'
  respiration <- x%v%'respiration'
  export <- x%v%'export'
                                        #recombine into the extended format
  import <- c(input,0,0,0)
  x <- cbind(flow,export,respiration,rep(0,nrow(flow)))
  x <- rbind(x,rep(0,length(import)),rep(0,length(import)),import)
                                        #make NA values zero if zero.na == TRUE
  if (zero.na){
    x[is.na(x)] = 0
  }
  return(x)
}
#' bal --- balances a flow model
#' INPUT = network model in extended format
#' OUTPUT = balanced model in extended format
#' NOTE: this is the work horse for balance.R
#' Original: M. Lau | July 2011
#' Re-written: M. Lau | 17Oct2013 
#' ---------------------------------------------------







#' bal --- balances a flow model INPUT = network model in extended format
#' OUTPUT = balanced model in extended format NOTE: this is the work horse for
#' balance.R Original: M. Lau | July 2011 Re-written: M. Lau | 17Oct2013
#' --------------------------------------------------- bal --- balances a flow
#' model INPUT = network model in extended format OUTPUT = balanced model in
#' extended format NOTE: this is the work horse for balance.R Original: M. Lau
#' | July 2011 Re-written: M. Lau | 17Oct2013
#' --------------------------------------------------- Subfunction for
#' Balancing by Either Inputs or Outputs
#' 
#' Dependency for the \code{balance} function.
#' 
#' 
#' @param T.star Extended, unbalanced matrix.
#' @param method Balance by inputs or outputs.
#' @return Returns an extended matrix for balancing by inputs or outputs.
#' @author Matthew K. Lau Stuart R. Borrett
#' @seealso \code{\link{balance}}
#' @references Fath, B.D. and S.R. Borrett. 2006. A MATLAB function for network
#' environ analysis. Environmental Modelling & Software 21:375-405.
bal <- function(T.star='matrix',method=c('input','output')){
  
  if (length(method > 1)){method <- method[1]}
  if (method == 'output'){T.star <- t(T.star)}  
                                        # transpose matrix to use output method
                                        #From Allesina and Bondavalli 2003
                                        #Step 1. Check balancing
                                        #Done in balance()
                                        #Step 2. Get F.star
  F.star <- T.star
  N <- nrow(T.star) - 3
  for (i in 1:N){
    for (j in 1:(N+3)){
      if (apply(T.star,1,sum)[i] == 0){
        F.star[i,j] <- 0
      }else{
        F.star[i,j] <- T.star[i,j] / apply(T.star,1,sum)[i]
      }
    }
  }
                                        #Step 3. Get R
  
  I.R <- diag(rep(1,N)) #F.star[1:N,1:N] identity matrix 
  R <- t(F.star[1:N,1:N]) - I.R
                                        #Step 4. Invert R
  R <- ginv(R)
                                        #Step 5. Multiply every rij by its corresponding input and change sign
  for (i in 1:nrow(R)){
    for (j in 1:nrow(R)){
      R[i,j] <- -(R[i,j] * (T.star[N+1,j]+T.star[N+2,j]+T.star[N+3,j]))
    }
  }
                                        #Step 6. Build U vector
  U <- apply(R,1,sum)
                                        #Step 7. Multiply F.star by corresponding U
  T.star.bal <- T.star
  for (i in 1:N){
    for (j in 1:(N+3)){
      T.star.bal[i,j] <- F.star[i,j] * U[i]
    }
  }
                                        #Final transposition if method == output
  if (method == 'output'){T.star.bal <- t(T.star.bal)}

  return(T.star.bal)

}
#' Balance Flow Network Models
#' 
#' Applies the methods of Allesina and Bondavalli (2003) for balancing flow
#' network models.
#' 
#' 
#' @param x A network object.
#' @param method Methods for model balancing, see Allesina and Bondavalli
#' (2003).
#' @param tol Percent error tolerance used in the steady state check prior to
#' balancing.
#' @return Returns a network object with a balanced flow network model.
#' @author Matthew K. Lau Stuart R. Borrett
#' @seealso \code{\link{bal}}
#' @references Allesina, S., Bondavalli, C., 2003. Steady state of ecosystem
#' flow networks: a comparison between balancing procedures. Ecological
#' Modelling 165(2-3):231-239.
#' @examples
#' 
#' 
#' 
#' data(troModels)
#' balance(troModels[[6]])
#' 
#' 
#' 
#' @export balance
balance <-
  function(x,method=c('AVG2','AVG','IO','OI','I','O'),tol=5){
                                        #Check for network class
  if (class(x) != 'network'){warning('x is not a network class object')}
  eT <- as.extended(x) #convert to extended format
  n <- network.size(x)
                                        #checks
  check <- ssCheck(x,tol)
  if (check){
    print('BALANCED',quote=FALSE);
    x%n%'balanced' = TRUE;
    return(x);
    stop
  }else{
    method <- method[1]
    print(method,quote= FALSE)
                                        #balancing
    if (method == 'AVG'){  ##Using the AVG method
      T.bal = 0.5 * (bal(eT,'input') + bal(eT,'output'))
    }else if (method == 'AVG2'){   ##Using the AVG2 method
      T.bal <- 0.5 *  (bal((0.5 * bal(eT,'output') + 0.5 * eT),'input')
                       + bal((0.5 * bal(eT,'input') + 0.5 * eT),'output'))
    }else if (method == 'IO'){   ##Using the IO method
      T.bal <- bal((0.5 * bal(eT,'input') + 0.5 * eT),'output')
    }else if (method == 'OI'){   ##Using the OI method
      T.bal <- bal((0.5 * bal(eT,'output') + 0.5 * eT),'input')
    }else if (method == 'I'){  # using the Input method
      T.bal <- bal(eT,'input')
    }else if (method == 'O'){
      T.bal <- bal(eT,'output')
    }else {warning('Unknown balancing method')}
                                        #convert balanced model into network class
    f <- T.bal[1:n,1:n] # create flow matrix
    set.edge.attribute(x,'flow',f[f>0])  # add flow weights to network edges
    x%v%'input' <- T.bal[(n+3),1:n]
    x%v%'export' <- T.bal[1:n,(n+1)]
    x%v%'respiration' <- T.bal[1:n,(n+2)]
    x%v%'output' <- (x%v%'export' + x%v%'respiration')
    x%v%'storage' <- x%v%'storage'
                                        #check for balancing and return output
    if (ssCheck(x)){
      x%n%'balanced' <- TRUE
      return(x)
    }else{
                                        #return false for unbalanced models
      warning('Model was not balanced.')}
      x%n%'balanced' <- FALSE
      return(x)
  }

}
#' Calculates the Ratio of Positive to Negative Elements in a Network
#' 
#' Dependent function for the enaUtility function.
#' 
#' 
#' @param x A matrix of flow values.
#' @return Returns the ratio of positive to negative elements in the flow
#' matrix.
#' @author Stuart R. Borrett
#' @seealso \code{\link{enaUtility}}
#' @references Fath, B.D. and S.R. Borrett. 2006. A MATLAB function for network
#' environ analysis. Environmental Modelling & Software 21:375-405.
bcratio <- function(x='matrix'){
	r=sum(x[x>0])/sum(abs(x[x<0]))  # creates the ratio of positive elements to negative elements.
  return(r)
}
#'## Cycle Analysis for Feeding Cycles
#'## Singh P.  | July 2014
#'## Algorithm Source : Ulanowicz 1991: A package for the Analysis of Ecosystem Flow Networks
#'## ---------------------------------------------

#' ## Cycle Analysis for Feeding Cycles ## Singh P.  | July 2014 ## Algorithm
#' Source : Ulanowicz 1991: A package for the Analysis of Ecosystem Flow
#' Networks ## --------------------------------------------- ## Cycle Analysis
#' for Feeding Cycles ## Singh P.  | July 2014 ## Algorithm Source : Ulanowicz
#' 1991: A package for the Analysis of Ecosystem Flow Networks ##
#' --------------------------------------------- Analysis of Feeding Cycles in
#' a Network
#' 
#' Performs the full cycle analysis on the living subset of the network based
#' on the algorithm described in Ulanowicz (1983) and implemented in NETWRK
#' 4.2b. It returns data.frames with details of the simple cycles and nexus,
#' vectors of Cycle distributions and Normalized distribution and matrices of
#' Residual Flows and Aggregated Cycles.
#' 
#' 
#' @param x a network object.  This includes all weighted flows into and out of
#' each node. It must also include the "Living" vector that identifies the
#' living (TRUE/FALSE) status of each node. Also, non-living nodes must be
#' placed at the end of the node vector. The function netOrder can be used to
#' reorder the network for this.
#' @return \item{Table.cycle}{data.frame that presents the details of the
#' simple cycles in the network. It contains "CYCLE" the cycle number, "NEXUS"
#' the nexus number corresponding to the cycle, "NODES" the nodes constituting
#' the cycle} \item{Table.nexus}{data.frame that presents the different nexuses
#' characterized by their corresponding weak arcs. It contains "NEXUS" the
#' nexus number, "CYCLES" the number of simple cycles present in that Nexus,
#' "W.arc.From" the starting node of the corresponding weak arc, "W.arc.To" the
#' ending node of the corresponding weak arc and "W.arc.Flow" the flow through
#' that weak arc} \item{CycleDist}{vector of the Cycle Distribution that gives
#' the flow which is cycling in loops of different sizes}
#' \item{NormDist}{vector of the Normalized Distribution i.e. the Cycle
#' Distribution normalized by the Total System Throughput for the network}
#' \item{ResidualFlows}{matrix of the straight-through (acyclic) flows in the
#' network} \item{AggregatedCycles}{matrix of the Aggregated Biogeochemical
#' Cycles in the network} \item{ns}{vector of the full cycle analysis based
#' network statistics. These include "NCYCS" the number of simple cycles
#' identified in the network, "NNEX" the number of the disjoint cycles or
#' number of Nexuses detected in the network and "CI" the cycling index of the
#' network.}
#' @note This function uses the same mechanism for analysis as used in the
#' enaCycle function but is restricted to the living nodes only.
#' 
#' Also, similar to the enaCycle function, if the number of cycles in a nexus
#' is more than 50, the "Table.cycle" has a blank line after 50 cycles followed
#' by the cycles for the next nexus.
#' 
#' The analysis requires all the non-living nodes to be placed at the end in
#' the network object.
#' @author Pawandeep Singh
#' @seealso \code{\link{enaTroAgg}, \link{enaCycle}, \link{netOrder}}
#' @references %% ~put references to the literature/web site here ~ Johnson,
#' D.B. 1975. Finding all the elementary circuits of a directed graph. SIAM J.
#' Comput. 4:77--84
#' 
#' Ulanowicz, R.E. 1983. Identifying the structure of cycling in ecosystems.
#' Methematical Biosciences 65:219--237
#' 
#' Ulanowicz, R.E. and Kay, J.J. 1991. A package for the analysis of ecosystem
#' flow networks. Environmental Software 6:131 -- 142.
cycliv <- function(x){

		 #Initials
    if(class(x)!='network') {stop("x is not a network class object")}
    web.all <- as.matrix(x,attrname="flow")
    y.all   <- x %v% "output"
    liv <- x %v% 'living'
    TPTS.all <- apply(web.all,1,sum)+y.all
    nl<-sum(liv)
    N<-length(liv)
    liv2<-rep(FALSE,N)
    liv2[1:nl]<-rep(TRUE,nl)
    if(identical(liv,liv2) == FALSE) {
    	stop('Non-living nodes must be at the end of the list.')
  }
    #-------------------------------
    N <- sum(liv)
    web <- web.all[1:N,1:N]
    y <- y.all[1:N]
    z<-(x %v% "input")[1:N]
    TPTS <- apply(web,2,sum)+z
    F <- web/TPTS

    TST <- sum(web)+sum(y)+sum(z)
    # -----------------------------
    ###-----------------------------------------------------------------
    df<-data.frame(NULL)
    df.cycle<-data.frame(NULL)
#'##-----------------------------------------------------------------

                                        #Zero Global Variables
    NFST <- NEXNUM <- NCYC <- 0
    CYCS <- rep(0,N)

#'##-----------------------------------------------------------------

                                        #Start primary repeat loop
    repeat {
                                        # Zero all local variables
        NNEX  <- 0
        TCYCS <- rep(0,N)
        TMP   <- web*0
                                        #-----------------------------------

                                        #Count cycle arcs and determine exit from return
        NFWD <- NULL
        NTEMP <-NULL
        for (ii in 1:N) {#do 200 ii=1,N
            NFWD <- rep(0,N)
            NFWD[ii] <- 1
            for (k in 1:(N-1)) {
                for (i in 1:N) {
                    if(NFWD[i]>0) {next}
                    for (j in 1:N) {
                        if (NFWD[j] <1) {next}
                        if (web[j,i]<=0) {next}
                        NFWD[i] <- 1
                        break
                    }
                }
            }
            NTEMP[ii] <- 0
            for (i in 1:N) {
                if ((NFWD[i]>0) && (web[i,ii]>0)) {NTEMP[ii] <- NTEMP[ii]+1}
            }
        }#200
                                        # NTEMP will give the no. of Cycles ending in each node
        NSTP <- 0
        MAP <- NULL
        for (i in 1:N) {
            NMAX <- -1
            for (j in 1:N) {
                if (NTEMP[j]<=NMAX) {next}
                NMAX <- NTEMP[j]
                JMAX <- j
            }
            if (NMAX > 0) {NSTP <- NSTP + 1}
            NTEMP[JMAX] <- -2
            MAP[i] <- JMAX
        }
        #print(c('NSTP',NSTP))


        ### Condition for breaking/Exiting the primary repeat loop
        ###----------------------------------------
        if (NSTP<=0) {break} #breaks the primary repeat loop
        ###----------------------------------------

                                        # Start the NSTP2 While loop for Critical Arc determination
        ###----------------------------------------------------------------------------------------------
        NSTP2 <- 0
        repeat {
            slf.loop<-FALSE
            ARCMIN <- 10^25 #arbitrary min arc
            for (ir in 1:NSTP) {
                for (ic in 1:NSTP) {
                    IRTP <- MAP[ir]  # IRTP is the node to be searched at the ir'th position
                    ICTP <- MAP[ic]  # ICTP ----------------------------------ic'th --------
                    if (web[IRTP,ICTP] <= 0) {next}
                    if (web[IRTP,ICTP] >= ARCMIN) {next}
                    ARCMIN <- web[IRTP,ICTP]
                    IMIN <- IRTP
                    IM <- ir
                    JMIN <- ICTP
                    JM <- ic
                }
            }
            #print(ARCMIN)
            #print(min(web[web>0]))
            ### Exit from while(NSTP2<=0) if slf.loop
            if (IMIN == JMIN) {
                slf.loop <- TRUE
                break    #-----------------------BREAK THE WHILE repeat LOOP
            }
                                        #Make sure at least one cycle contains the current smallest arc web[IMIN,JMIN]
            NHALF <- (N/2)+1
            NFWD <- rep(0,N)
            NFWD[JMIN] <- 1
            ### find nodes from JMIN in fwd dirctn
            for (k in 1:NHALF) {
                for (i in 1:N) {
                    if (NFWD[i] >0) {next}
                    for (j in 1:N) {
                        if (NFWD[j] < 1) {next}
                        if (web[j,i] <= 0) {next}
                        NFWD[i] <- 1
                        break
                    }
                }
            }
            ### find nodes from IMIN in bckwd dirctn
            NODE <- rep(0,N)
            NODE[IMIN] <- 1
            for (k in 1:NHALF) {
                for (i in 1:N) {
                    if (NODE[i] > 0) {next}
                    for (j in 1:N) {
                        if (NODE[j] < 1) {next}
                        if (web[i,j] <= 0) {next}
                        NODE[i] <- 1
                        break
                    }
                }
            }
            ### find Common nodes aka members of the NEXUS
            NSTP2 <- 0  ###### ?????
            MAP2 <- rep(0,NSTP)
            for (i in 1:NSTP) {
                if ((NFWD[MAP[i]] <= 0) | (NODE[MAP[i]] <= 0)) {next}
                NSTP2 <- NSTP2 + 1
                MAP2[NSTP2] <- MAP[i]
            }
                                        # reorder mapping for IMIN and JMIN to come 1st and 2nd
            NFWD <- MAP2
            MAP2[1] <- IMIN
            MAP2[2] <- JMIN ##### SHORTER WAY POSSIBLE ####
            if (NSTP2 > 2) {
                INDX <- 2
                for (i in 1:NSTP2) {
                    if ((NFWD[i] == IMIN)||(NFWD[i]==JMIN)) {next}
                    INDX <- INDX+1
                    MAP2[INDX] <- NFWD[i]
                }
            }
            if (NSTP2 > 0) {break}
            web[IMIN,JMIN] = -web[IMIN,JMIN]
        } ###End of While Repeat NSTP2<=0
        ###----------------------------------------------------------------------------------------

                                        #IF the Critical arc is self loop
                                       #---------------------------------
        if(slf.loop == TRUE) {
            slf.loop <- FALSE
        	NCYC <- NCYC+1
            NNEX <- NNEX+1
            CYCS[1]<- CYCS[1]+web[IMIN,JMIN]
            WKARC<-F[IMIN,JMIN]*TPTS[IMIN]
            curr.slf.cyc<-noquote(c(NCYC,'.','(',IMIN,JMIN,')'))
            #print(curr.slf.cyc)
            this.cycle <- rep(NA,N)
            this.cycle[1:2] <- c(IMIN,JMIN)
            newcycle<-c(NCYC,(NEXNUM+1),this.cycle)
            df.cycle<-rbind(df.cycle,newcycle)

            web[IMIN,JMIN] <- 0
            NEXNUM <- NEXNUM+1
            curr.nexus <- c(NEXNUM,NNEX,IMIN,JMIN,WKARC)

            df<-rbind(df,curr.nexus) ### ('NEXUS', 'Cycles', weakarc, fromnode, tonode) ####################### df

            NFST <- 1
        }#End of if(slf.loop==TRUE)#

                                        #Begin Search for NEXUS defined by web[IMIN,JMIN] if not a self loop
        ###----------------------------------------------------------------------------------------

        else {
            WHOLE <- 0
                                        # Backtrack Routine Starts --------------
             # Backtrack Routine Starts --------------
            ### Initialize Node and Level
            LEVEL  <- 2
            NODE[1]<- 1
            NODE[2]<- 2
            skip.con.adv <- FALSE
            nex.com <- FALSE

            ### 2 Repeats start. rep1,2.
            repeat { #rep1
                repeat { #rep2
                    skip.con.chk <- FALSE

             ## Adv to nxt levels
                    if(skip.con.adv==FALSE) {
                        LM1         <- LEVEL
                        LEVEL       <- LEVEL+1
                        NODE[LEVEL] <- 1
                    }
#'## 2 Repeats start. rep3,4
                    repeat { #rep3
                        repeat { #rep4
                            ## Check for conn. b/w nodes at prsnt levels
                            if(skip.con.adv==FALSE){
                                NZ1  <- NODE[LM1]
                                KROW <- MAP2[NZ1]
                                NZ2  <- NODE[LEVEL]
                                KCOL <- MAP2[NZ2]
                                conn.chk <- FALSE
                                ## break rep4 and rep3 if conn exists
                                if(web[KROW,KCOL]>0 && skip.con.chk==FALSE) {
                                    conn.chk <- TRUE
                                    break #brk rep4
                                }
                            }
                            ## try next node in nxt level
                            NODE[LEVEL] <- NODE[LEVEL]+1
                            skip.con.adv <- FALSE
                            skip.con.chk <- FALSE
                            ## break rep4 after all levels are checked(NODE[level]>NSTP2)
                            if(NODE[LEVEL]>NSTP2) {break} #brk rep4
                        }#end of rep4
                        if(conn.chk==TRUE) {
                            conn.chk <- FALSE
                            break #rep3
                        }
#'## Backtrack to prev. level
                        LEVEL <- LEVEL-1
                        LM1   <- LEVEL-1
                        ## if further backtracking is impossible,
                                        #end search under weak arc. nexus complete
                        nex.com<-FALSE
                        if(LEVEL <= 2){
                            nex.com <- TRUE
                            break #break rep3
                        }
                        else {skip.con.chk<-TRUE} #goes to #420 to inc NODE[LEVEL]
                    }#end of rep3
                    if(nex.com==TRUE) {break} #break rep2
                    ##break rep2 if this conn completes cycle
                    if(NODE[LEVEL]==1) {break} #brk rep2
                    skip.con.adv <- FALSE
                    for (k in 1:LM1) {if (NODE[LEVEL] == NODE[k]) {skip.con.adv <- TRUE}}

                }#end of rep2
                if(nex.com==TRUE) {break} #brk rep1
                                        # -----------------BR Ends --------------
                                        # Calculate circuit prob
                WEIGHT <- 1
                for (kk in 1:LM1) {
                    KKP1 <- kk+1
                    KROW <- NODE[kk]
                    KCOL <- NODE[KKP1]
                    KROW <- MAP2[KROW]
                    KCOL <- MAP2[KCOL]
                    WEIGHT <- WEIGHT*F[KROW,KCOL]
                }
                                        # Add this weight
                for (kk in 1:LM1) {
                    KKP1 <- kk+1
                    KROW <- NODE[kk]
                    KCOL <- NODE[KKP1]
                    KROW <- MAP2[KROW]
                    KCOL <- MAP2[KCOL]
                    # Also, add this amount to the cycle distributions
                    TCYCS[LM1] <- TCYCS[LM1]+WEIGHT
                    TMP[KROW,KCOL] <- TMP[KROW,KCOL]+WEIGHT
                }
                                        # Report this cycle
                NNEX <- NNEX+1
                KTRY <- NNEX%%5000
                curr.prog <- noquote(c(NNEX,'Nexus cycles and Counting'))
                if(KTRY==0) {print(curr.prog)}
                NCYC <- NCYC+1
                if(NNEX>50){
                	skip.con.adv <- TRUE
                	next
                }
                L0 <- LM1+1
                for (kk in 1:L0) {
                    NTMP <- NODE[kk]
                    NTEMP[kk] <- MAP2[NTMP]
                }
                #curr.cycle <- noquote(c(NCYC,'.',NTEMP[1:L0]))
                #print(curr.cycle)
                this.cycle <- NTEMP
                this.cycle[(L0+1):N]<-NA
                newcycle<-c(NCYC,(NEXNUM+1),this.cycle)
                df.cycle<-rbind(df.cycle,newcycle)            #################------------------------------------df.cycle
                if(NNEX==50) {df.cycle<-rbind(df.cycle,rep(NA,N+2))}
                skip.con.adv<-TRUE
            }#end of rep1
                                        #-----------------NEXUS COMPLETED---NEXUS REPEAT(rep1) ENDS HERE



                                        # Report this NEXUS
            WKARC=F[IMIN,JMIN]*TPTS[IMIN]
            NEXNUM=NEXNUM+1
            curr.nexus <- c(NEXNUM,NNEX,IMIN,JMIN,WKARC)
            df<-rbind(df,curr.nexus) ### ('NEXUS', 'Cycles', weakarc, fromnode, tonode) ####################### df
            #print(curr.nexus)
            # -------------------------------------
            # Normalize Probability Matrix & Subtract proper amounts from web
            PIVOT <- TMP[IMIN,JMIN]
            if(PIVOT<=0){print('Error in Normalizing Nexus Weights')}
            for(i in 1:N) {
                for(j in 1:N) {
                    if(web[i,j]<=0) {next}
                    web[i,j] <- web[i,j]-((TMP[i,j]/PIVOT)*ARCMIN)
                }
            }

            # Add proper amts to cycle distributions
            for(i in 1:N){CYCS[i]<-CYCS[i]+((TCYCS[i]/PIVOT)*ARCMIN)}               ############################## ERR.CHK CYCS depends on TCYCS, PIVOT, ARCMIN

            # Zero weak arc
            web[IMIN,JMIN] <- 0
            NFST <- 1
            if(WHOLE>1.00001) {print(c('Bad Sum Check = ',WHOLE))}


        }#End of else i.e. if not a slf.loop
    }#End of primary repeat loop

                                              # Report the Overall Results

    ### FIRST, UNCOVER ANY LINKS "HIDDEN" DURING SEARCH.
    web=abs(web)
    if(NFST!=0) {
        #cyc.rem <- noquote((c('A total of ',NCYC,'Cycles removed')))
        #print(cyc.rem)
        #print('Cycle Distributions')
        #print(CYCS)
        cycs<-CYCS
        CYC  <- sum(CYCS)
        CYCS <- (CYCS/TST)
        #print('Normalized Distribution')
        #print(CYCS)
        TEMP <- CYC/TST
        #print(c('cycling index is',TEMP))
        ResidualFlows<-web
        AggregatedCycles<-(as.matrix(x, attrname = "flow")[1:N,1:N]) - ResidualFlows
        colnames(df)<-c('NEXUS', 'Cycles','From','To', 'Weak_arc')
        colnames(df.cycle)<-rep(' ',(N+2))
        colnames(df.cycle)[1:3]<-c('CYCLE','NEXUS','NODES')
        df.cycle[is.na(df.cycle)==TRUE]<- ' '
        NCYCS<-NCYC; NNEX<-NEXNUM; CI<-TEMP
        ns <- cbind(NCYCS, NNEX, CI)
        out <- list(Table.cycle=df.cycle,Table.nexus=df, CycleDist = cycs, NormDist=CYCS, ResidualFlows=web, AggregatedCycles=AggregatedCycles, ns=ns)
        return(out)
    }#end of if (NFST!=0)
    else {
        NCYCS<-NCYC; NNEX<-NEXNUM; CI <- 0
        ns <- cbind(NCYCS,NNEX,CI)
        out <- list(ResidualFlows=web, ns=ns)
        return(out)
      }


}#END OF FUNCTION
#' Calculates the Eigen Centrality of a Network
#' 
#' Calculates the eigen centrality of a network.
#' 
#' 
#' @param x A matrix defining a network graph.
#' @return Returns the eigen based centrality of the network.
#' @author Stuart R. Borrett Matthew K. Lau
#' @references Bonacich, P., 1987. Power and centrality: a family of measures.
#' American Journal of Sociology 92: 1170-1182.
#' @export eigenCentrality
eigenCentrality <- function(x='matrix'){
  if (class(x) != 'matrix'){warning('x is not a matrix class object')}
                                        # find dominant eigenvector of x
  EVCin <- abs(eigen(x)$vectors[,1])
  EVCin <- EVCin/sum(EVCin)           # normalize by sum
                                        # find dominant eigenvector of x transpose
  EVCout <- abs(eigen(t(x))$vectors[,1])  
  EVCout <- EVCout/sum(EVCout)        # normalize by sum  
  AEVC <- (EVCin + EVCout)/2          # find average eigenvector centrality
  
  return(list('EVCin'=EVCin,'EVCout'=EVCout,'AEVC'=AEVC))
}
#' enaAll --- Conduct all ecological network analyses
#' INPUT = network object
#' OUTPUT = list of analytical output
#' 
#' M. Lau | May 2013
#' ---------------------------------------------------







#' enaAll --- Conduct all ecological network analyses INPUT = network object
#' OUTPUT = list of analytical output
#' 
#' M. Lau | May 2013 --------------------------------------------------- enaAll
#' --- Conduct all ecological network analyses INPUT = network object OUTPUT =
#' list of analytical output
#' 
#' M. Lau | May 2013 ---------------------------------------------------
#' Conduct All Major ENA
#' 
#' Conducts all major ENA with default settings and returns the output as a
#' named list.
#' 
#' @param x A network object.
#' @return \item{ascendency}{enaAscendency} \item{control}{enaControl}
#' \item{environ}{enaEnviron} \item{flow}{enaFlow} \item{mti}{enaMTI}
#' \item{storage}{enaStorage} \item{structure}{enaStructure}
#' \item{utility}{enaUtility with eigen.check=FALSE}
#' @author Matthew K. Lau Stuart R. Borrett
#' @seealso
#' \code{\link{enaAscendency},\link{enaControl},\link{enaEnviron},\link{enaFlow},\link{enaMTI},\link{enaStorage},\link{enaUtility}}
#' @examples
#' 
#' 
#' 
#' data(troModels)
#' output = enaAll(troModels[[6]])
#' names(output)
#' 
#' 
#' 
#' @export enaAll
enaAll <- function(x = 'network object'){
  out <- list(ascendency = enaAscendency(x),
              control = enaControl(x),
              environ = enaEnviron(x),
              flow = enaFlow(x),
              mti = enaMTI(x),
              storage = enaStorage(x),
              structure = enaStructure(x),
              utility = enaUtility(x,eigen.check=FALSE))
  return(out)
}
#' enaAscendency --- calculates the ascendency statistics
#' of Ulanowicz
#' INPUT = network object
#' OUTPUT = matrix of ascendency statistics
#'
#' D. Hines | December 2011
#' S.R. Borrett | May 2016 - updates
#' ---------------------------------------------------







#' enaAscendency --- calculates the ascendency statistics of Ulanowicz INPUT =
#' network object OUTPUT = matrix of ascendency statistics
#' 
#' D. Hines | December 2011 S.R. Borrett | May 2016 - updates
#' --------------------------------------------------- enaAscendency ---
#' calculates the ascendency statistics of Ulanowicz INPUT = network object
#' OUTPUT = matrix of ascendency statistics
#' 
#' D. Hines | December 2011 S.R. Borrett | May 2016 - updates
#' --------------------------------------------------- Calculates the
#' Ascendency of an Ecological Network
#' 
#' Calculates the average mutual information (AMI), ascendency, overhead, and
#' capacity of input-output networks.  It also returns the ratios of ascendency
#' and overhead to capacity. These metrics describe the organization of flow in
#' an ecological network (Ulanowicz 1997).
#' 
#' @param x A network object.
#' @return \item{H}{Total flow diversity (MacArthur 1955).  Uses the Shannon
#' Information measure (aka Boltzmann entropy) applied to the individual flows.
#' } \item{AMI}{Returns the Average Mutual Information (AMI) in a network. AMI
#' provides a measure of the constraints placed on a given peice of energy
#' matter moving through a network (Patricio et al. 2006) } \item{Hr}{Residual
#' uncertainty that remains about the flow distribution once the ecosystem
#' structure is specified (Hr = H - AMI). } \item{ASC}{Returns the ascendnecy
#' of a network.  Ascendency is a scaled form of AMI relative to the total
#' system throughput (Ulanowicz 1997; 2004).  Total system throughput is the
#' sum of all activity in a network (Kay et al. 1989).} \item{OH}{Returns the
#' overhead of a network.  Overhead is the proportion of the capacity in a
#' network that is not used as ascendency (Ulanowicz 2004).} \item{CAP}{Returns
#' the capacity of a network.  Capacity is defined as the sum of ascendency and
#' overhead (Ulanowicz 2004).} \item{ACS.CAP}{Returns the proportion of
#' capacity used by ascendency.} \item{OH.CAP}{Returns the proportion of
#' capacity used by overhead.} \item{robustness}{Returns the robustness of the
#' network.} \item{ELD}{Returns the Effective Link Density of the network(c)
#' (Ulanowicz et al. 2014).} \item{TD}{Returns the Trophic Depth of the
#' network(r) (Ulanowicz et al. 2014).} \item{A.input}{Returns the input
#' ascendnecy of a network.} \item{A.internal}{Returns the internal ascendnecy
#' of a network.} \item{A.export}{Returns the export ascendnecy of a network.}
#' \item{A.respiration}{Returns the respiration ascendnecy of a network.}
#' \item{OH.input}{Returns the input overhead of a network.}
#' \item{OH.internal}{Returns the internal overhead of a network.}
#' \item{OH.export}{Returns the export overhead of a network.}
#' \item{OH.respiration}{Returns the respiration overhead of a network.}
#' \item{CAP.input}{Returns the input capacity of a network.}
#' \item{CAP.internal}{Returns the internal capacity of a network.}
#' \item{CAP.export}{Returns the export capacity of a network.}
#' \item{CAP.respiration}{Returns the respiration capacity of a network.}
#' @note This and other Ulanowicz school functions require that export and
#' respiration components of output be separately quantified.
#' @author David E. Hines Matthew K. Lau Stuart R. Borrett
#' @seealso
#' \code{\link{read.scor},\link{read.wand},\link{enaStorage},\link{enaUtility}}
#' @references Kay, J.J., Graham, L.A., Ulanowicz, R.E., 1989. A detailed guide
#' to network analysis. p. 15-61 In: Wulff, F., Field, J.G., Man, K.H. (eds.)
#' Network analysis in marine ecology. Coastal Estuarine Study Serries.
#' Springer-Verlag, Berlin.
#' 
#' Patrico, J., Ulanowicz, R.E., Pardal, M.A., Marques J.C., 2004. Ascendency
#' as an ecological indicator: a case study of estuarine pulse eutrophication.
#' Estuar. Coast Shelf S. 60, 23-35.
#' 
#' Ulanowicz, R.E. and Norden, J.S., 1990. Symmetrical overhead in flow
#' networks. International Journal of Systems Science, 21(2), pp.429-437.
#' 
#' Ulanowicz, R.E., 1997. Ecology, The Ascendent Perspective. Columbia
#' University Press, New York.
#' 
#' Ulanowicz, R.E., 2004. Quantitative methods for ecological network analysis.
#' Comput. Biol. Chem. 28, 321-33
#' 
#' Ulanowicz, R.E., Holt, R.D., Barfield, M., 2014. Limits on ecosystem trophic
#' complexity: insights from ecological network analysis. Ecology Letters
#' 17:127-136
#' @export enaAscendency
#' @examples
#' data(troModels)
#' enaAscendency(troModels[[6]])
enaAscendency <- function(x='network object'){
    if (class(x) != 'network'){warning('x is not a network class object')}


    if (any(is.na(x%v%'export'))){
        warning('Export data is absent from the model.')
    }
    if(any(is.na(x%v%'respiration'))){
           warning('Respiration data is absent from the model.')
   }

#'####### set initial conditions for calculations #########
  T.ulan <- as.extended(x)
  N <- ncol(T.ulan) # set up N
  r.td <- c.ld <- t.ulan <- ami <- mat.or.vec(N,N) # initialize ascendency matrix
  oh <- mat.or.vec(N,N) # initialize overhead matrix
  cap <- mat.or.vec(N,N) # initialize capacity matrix
                                        #calculate total system throughPUT
  TSTp <- sum(T.ulan)

#'## calculate H & CAPACITY  #######################################
  #' H = Total Flow Diversity

  h <- T.ulan/sum(T.ulan)
  h2 <- log2(h)
  h2[!is.finite(h2)] <- 0
  H = - sum(h * h2)   # Total Flow Diversity

  CAP <- H * TSTp     # Capactity

#'################### calculate AMI  #######################
  #' AMI = Average Mutual Informaiton

                                        # loop through T.ulan to calculate AMI
  for (i in 1:N){
    for (j in 1:N){
      if (T.ulan[i,j] == 0){
        ami[i,j] <- 0
      }else{
        ami[i,j] <- T.ulan[i,j]/TSTp * log2((T.ulan[i,j]*TSTp)/(sum(T.ulan[,j])*sum(T.ulan[i,])))
      }
    }
  }

  AMI <- sum(ami)

#'################ calculate ascendency ###################

  ASC <- TSTp * AMI

#'################ calculate residual diversity  ####################

  Hr <- H - AMI

#'################ calculate overhead  ####################

                                        # loop through T.ulan to calculate overhead
  for (i in 1:N){
    for (j in 1:N){
      if (T.ulan[i,j] == 0){
        oh[i,j] <- 0
      }else{
        oh[i,j] <- T.ulan[i,j] * log2((T.ulan[i,j]^2)/(sum(T.ulan[,j])*sum(T.ulan[i,])))
      }
    }
  }

  OH <- -sum(oh)


#'################### calculate ratios ####################

                                        # ratio for ascendency/capacity
  ASC.CAP <- ASC/CAP

                                        #ratio for overhead/capacity
  OH.CAP <- OH/CAP

                                        #confirm ratios sum to 1
  ASC.OH.RSUM <- ASC.CAP + OH.CAP

  robustness = -1 * ASC.CAP * log(ASC.CAP)  # robustness from Ulanowicz 2009; Fath 2014



  ################# Calculating Effective Link Density and Trophic Depth ########
  ## Calculate t.ulan 't'

  for (i in 1:N) {
        for (j in 1:N) {
            if (T.ulan[i, j] == 0) {
                t.ulan[i, j] <- 0
            }
            else {
                t.ulan[i, j] <- T.ulan[i,j]/TSTp
            }
        }
    }

    ## Effective Link Density (c)
    for (i in 1:N) {
        for (j in 1:N) {
            if (t.ulan[i, j] == 0) {
                c.ld[i, j] <- 1
            }
            else {
                c.ld[i, j] <- (sqrt(sum(t.ulan[i,])*sum(t.ulan[,j]))/t.ulan[i,j])^(t.ulan[i,j])
            }
        }
    }
    C.LD <- prod(c.ld)

    ## Trophic Depth (r)
    for (i in 1:N) {
        for (j in 1:N) {
            if (t.ulan[i, j] == 0) {
                r.td[i, j] <- 1
            }
            else {
                r.td[i, j] <- (t.ulan[i,j]/(sum(t.ulan[i,])*sum(t.ulan[,j])))^(t.ulan[i,j])
            }
        }
    }
    R.TD <- prod(r.td)

    ELD <- C.LD
    TD <- R.TD
    ##############################################################################


  #'#####################################################################
  # tetrad partition of A, C, O -> input, internal, export, respiration
  #'#####################################################################
  n <- N - 3

  # -- ASCENDENCY --
  # input
  tmp = 0
  for(j in 1:n){
      tmp[j] <-  T.ulan[(n+3),j]  * log2(  (T.ulan[(n+3),j] * TSTp)/( sum(T.ulan[(n+3),]) * sum(T.ulan[,j]) ) )
  }
  A.input <- sum(tmp[!is.nan(tmp)])   # have to remove NaN values

  # internal
  tmp=mat.or.vec(n,n)
  for(i in 1:n){
      for(j in 1:n){
          tmp[i,j]= T.ulan[i,j] * log2(  (T.ulan[i,j] * TSTp) /( sum(T.ulan[i,]) * sum( T.ulan[,j]) ))
      }
  }
  A.internal <- sum(tmp[!is.nan(tmp)])

  # Exports
  tmp <- 0
  for(i in 1:n){
      tmp[i] = T.ulan[i,(n+1)] * log2( (T.ulan[i,(n+1)] * TSTp) / ( sum(T.ulan[,(n+1)]) * sum(T.ulan[i,]) ) )
  }
  A.export <- sum(tmp[!is.nan(tmp)])

  # Respiration
  tmp <- 0
  for(i in 1:n){
      tmp[i] = T.ulan[i,(n+2)] * log2( (T.ulan[i,(n+2)] * TSTp) / ( sum(T.ulan[,(n+2)]) * sum(T.ulan[i,]) ) )
  }
  A.respiration <- sum(tmp[!is.nan(tmp)])

  # -- OVERHEAD --
  # input
  tmp = 0
  for(j in 1:n){
      tmp[j] <-  -1 * T.ulan[(n+3),j]  * log2(  (T.ulan[(n+3),j]^2 )/( sum(T.ulan[(n+3),]) * sum(T.ulan[,j]) ) )
  }
  OH.input <- sum(tmp[!is.nan(tmp)])   # have to remove NaN values

  # internal
  tmp=mat.or.vec(n,n)
  for(i in 1:n){
      for(j in 1:n){
          tmp[i,j]= -1 * T.ulan[i,j] * log2(  (T.ulan[i,j]^2) /( sum(T.ulan[i,]) * sum( T.ulan[,j]) ))
      }
  }
  OH.internal <- sum(tmp[!is.nan(tmp)])

  # export
  tmp <- 0
  for(i in 1:n){
      tmp[i] = -1 * T.ulan[i,(n+1)] * log2( (T.ulan[i,(n+1)]^2) / ( sum(T.ulan[,(n+1)]) * sum(T.ulan[i,]) ) )
  }
  OH.export <- sum(tmp[!is.nan(tmp)])

  # respriation
  tmp <- 0
  for(i in 1:n){
      tmp[i] = -1 *  T.ulan[i,(n+2)] * log2( (T.ulan[i,(n+2)]^2) / ( sum(T.ulan[,(n+2)]) * sum(T.ulan[i,]) ) )
  }
  OH.respiration <- sum(tmp[!is.nan(tmp)])

  # -- CAPACITY --
  CAP.input = A.input + OH.input
  CAP.internal = A.internal + OH.internal
  CAP.export = A.export + OH.export
  CAP.respiration = A.respiration + OH.respiration


  ns <- cbind(H, AMI, Hr, CAP, ASC, OH, ASC.CAP, OH.CAP,
              robustness, ELD, TD,
              A.input, A.internal, A.export, A.respiration,
              OH.input, OH.internal, OH.export, OH.respiration,
              CAP.input, CAP.internal, CAP.export, CAP.respiration)
                                        #
  return(ns)

}
#' enaControl --- control analyses
#' INPUT = network object
#' OUTPUT = list of control statistics
#' M. Lau | July 2011
#' P. Singh | Update Summer 2013
#' S.R. Borrett | Update March 2016

#' ---------------------------------------------------







#' enaControl --- control analyses INPUT = network object OUTPUT = list of
#' control statistics M. Lau | July 2011 P. Singh | Update Summer 2013 S.R.
#' Borrett | Update March 2016
#' --------------------------------------------------- enaControl --- control
#' analyses INPUT = network object OUTPUT = list of control statistics M. Lau |
#' July 2011 P. Singh | Update Summer 2013 S.R. Borrett | Update March 2016
#' --------------------------------------------------- Control Analyses of
#' Ecological Networks
#' 
#' Analyses for analyzing the control amongst the nodes in ecological networks.
#' 
#' 
#' @param x A network object.
#' @param zero.na Makes undefined (NA) values zero.
#' @param balance.override Turns off balancing and checks of network balance.
#' @return \item{CN}{Control matrix using flow values.} \item{CQ}{Control
#' matrix using storage values.} \item{CR}{Schramski Control Ratio Matrix}
#' \item{CD}{Schramski Control Difference Matrix} \item{CA}{Control Allocation
#' Matrix} \item{CDep}{Control Dependency Matrix} \item{sc}{Schramski System
#' Control vector} \item{scp}{Schramski system control vector as percent of
#' total control} \item{ns}{vector of network-level summary statistics}
#' @author Matthew K. Lau Stuart R. Borrett Pawandeep Singh
#' @seealso \code{\link{enaStorage}}
#' @references Fath, B. D., Borrett, S. R. 2006. A MATLAB function for Network
#' Environ Analysis.  Environmental Modelling & Software 21:375-405
#' 
#' Schramski, J.R., Gattie, D.K., Patten, B.C., Borrett S.R., Fath, B.D.,
#' Thomas, C.R., and Whipple, S.J. 2006. Indirect effects and distributed
#' control in ecosystems: Distributed control in the environ networks of a
#' seven compartment model of nitrogen flow in the Neuse River Estuary, USA
#' Steady-state analysis. Ecological Modelling 194:189-201
#' 
#' Schramski, J.R., Gattie, D.K., Patten, B.C., Borrett S.R., Fath, B.D., and
#' Whipple, S.J. 2007. Indirect effects and distributed control in ecosystems:
#' Distributed control in the environ networks of a seven compartment model of
#' nitrogen flow in the Neuse River Estuary, USA Time series analysis.
#' Ecological Modelling 206:18-30
#' 
#' Chen, S., Fath, B.D., Chen, B. 2011. Information-based network environ
#' analysis: a system perspective for ecologcial risk assessment.  Ecol. Ind.
#' 11:1664-1672.
#' 
#' Chen, S. and Chen, B. 2015. Urban energy consumption: Different insights
#' from energy flow analysis, input-output analysis and ecological network
#' analysis.  Applied Energy 138:99-107.
#' @examples
#' 
#' 
#' 
#' data(troModels)
#' enaControl(troModels[[6]])
#' 
#' 
#' 
#' @export enaControl
enaControl <- function(x, zero.na=TRUE,balance.override=FALSE){
                                        #Check for network class
    if (class(x) != 'network'){warning('x is not a network class object')}
                                        #Check for balancing
    if (balance.override){}else{
        if (any(list.network.attributes(x) == 'balanced') == FALSE){x%n%'balanced' <- ssCheck(x)}
        if (x%n%'balanced' == FALSE){warning('Model is not balanced'); stop}
    }

    u <- unpack(x)  # unpack the model
    F <- enaFlow(x, balance.override=balance.override) # perform flow analysis

    # Calculate the control matrix - FLOW
    CN <- F$N / t(F$NP)
    i <- which(CN<1,arr.ind=TRUE)
    j <- which(!(CN<1),arr.ind=TRUE)
    CN[i] <- 1 - CN[i]
    CN[j] <- 0
    if (zero.na){
        CN[!is.finite(CN)] <- 0
    }

    # Calculate the control matrix - STORAGE
    # H:  Storage and Flow provide the same answer
    if(!any(is.na(x%v%'storage'))){
        S <- enaStorage(x,balance.override=balance.override)
        CQ <- S$Q / t(S$QP)
        i <- which(CQ<1,arr.ind=TRUE)
        j <- which(!(CQ<1),arr.ind=TRUE)
        CQ[i] <- 1 - CQ[i]
        CQ[j] <- 0
        if (zero.na){
            CQ[!is.finite(CQ)] <- 0
        }
    } else {
        CQ <- NA
    }

        # Schramski Control Measures (2006, 2007)

    eta <- t(t(F$N) / F$T)  # row-to-column
    CD <- eta - t(eta)      # control difference
    CR <- CD/pmax(eta,t(eta))  # control ratio
    sc <- apply(CD,1,sum)  # system control vector
    psc <- sc/(sum(abs(sc)/2)) * 100 # percent system control vector

    TSC <- sum(abs(sc)/2)
    ns <- c("TSC"=TSC)

    # Control Allocation and Control Dependence
    # Chen et al. 2011; Chen and Chen 2015

    d <- F$N - t(F$NP)  # difference
    d[d<0] = 0  # remove negative values
    cd.r <- apply(d,1,sum)
    cd.c <- apply(d,2,sum)

    CA <- ginv(diag(cd.r)) %*% d   # control allocation matrix
    CDep <- ginv(diag(cd.c)) %*% d   # control depedency matrix


    orient <- get.orient()
    if (orient == 'school'){
        CN <- t(CN)
        CQ <- t(CQ)
        CR <- t(CR)
        CD <- t(CD)
        CA <- t(CA)
        CDep <- t(CDep)

    }

    return(list("CN"=CN,"CQ"=CQ,"CD"=CD,"CR"=CR, "CA"=CA, "CDep"=CDep,
                "sc"=sc,"psc"=psc, "ns"=ns))
}

#'## NETWRK's Full Cycle Analysis
#'## Singh P. | July 2014
#'## Algorithm Source : Ulanowicz 1991: A package for the Analysis of Ecosystem Flow Networks
#'## -----------------------------------------------








#' ## NETWRK's Full Cycle Analysis ## Singh P. | July 2014 ## Algorithm Source
#' : Ulanowicz 1991: A package for the Analysis of Ecosystem Flow Networks ##
#' ----------------------------------------------- ## NETWRK's Full Cycle
#' Analysis ## Singh P. | July 2014 ## Algorithm Source : Ulanowicz 1991: A
#' package for the Analysis of Ecosystem Flow Networks ##
#' ----------------------------------------------- Full Cycle Analysis of
#' Ecological Networks
#' 
#' It performs the full cycle analysis on the network based on the algorithm
#' described in Ulanowicz (1983) and implemented in NETWRK 4.2b. It returns
#' data.frames with details of the simple cycles and nexus, vectors of Cycle
#' distributions and Normalized distribution and matrices of Residual Flows and
#' Aggregated Cycles.
#' 
#' 
#' @param x a network object.  This includes all weighted flows into and out of
#' each node.
#' @return \item{Table.cycle}{data.frame that presents the details of the
#' simple cycles in the network. It contains "CYCLE" the cycle number, "NEXUS"
#' the nexus number corresponding to the cycle, "NODES" the nodes constituting
#' the cycle} \item{Table.nexus}{data.frame that presents the different nexuses
#' characterized by their corresponding weak arcs. It contains "NEXUS" the
#' nexus number, "CYCLES" the number of simple cycles present in that Nexus,
#' "W.arc.From" the starting node of the corresponding weak arc, "W.arc.To" the
#' ending node of the corresponding weak arc and "W.arc.Flow" the flow through
#' that weak arc} \item{CycleDist}{vector of the Cycle Distribution that gives
#' the flow cycling in loops of different sizes} \item{NormDist}{vector of the
#' Normalized Distribution i.e. the Cycle Distribution normalized by the Total
#' System Throughput of the system} \item{ResidualFlows}{matrix of the
#' straight-through (acyclic) flows in the network}
#' \item{AggregatedCycles}{matrix of the Aggregated Biogeochemical Cycles in
#' the network} \item{ns}{vector of the full cycle analysis based network
#' statistics. These include "NCYCS" the number of simple cycles identified in
#' the network, "NNEX" the number of the disjoint cycles or number of Nexuses
#' detected in the network and "CI" the cycling index of the network.}
#' @note The "NODES" in "Table.cycle" are arranged such that the weak arc for
#' the nexus is the arc between the first two nodes of the cycle.
#' 
#' This function uses the backtracking procedure for the identification of
#' simple cycles, which are cycles that cross a node only once except the
#' starting node. The backtracking process is a depth-first search algorithm.
#' 
#' In the data.frame "Table.cycle", if the number of cycles in a nexus is more
#' than 50, then a blank line is displayed after 50 cycles of the nexus,
#' followed by the cycles of the next nexus.
#' 
#' The results of the analysis of Feeding Cycles can be obtained as a byproduct
#' of the enaTroAgg function that analyzes the trophic dynamics of a network.
#' 
#' At every multiple of 5000 cycles in a nexus, the program prints an
#' indication for the user to know that it is still running.
#' @author Pawandeep Singh
#' @seealso \code{\link{enaTroAgg}}
#' @references %% ~put references to the literature/web site here ~ Johnson,
#' D.B. 1975. Finding all the elementary circuits of a directed graph. SIAM J.
#' Comput. 4:77--84
#' 
#' Ulanowicz, R.E. 1983. Identifying the structure of cycling in ecosystems.
#' Methematical Biosciences 65:219--237
#' 
#' Ulanowicz, R.E. and Kay, J.J. 1991. A package for the analysis of ecosystem
#' flow networks. Environmental Software 6:131 -- 142.
#' @export enaCycle
#' @examples
#' 
#' 
#' 
#' data(troModels)
#' cyc6 <- enaCycle(troModels[[6]])
#' attributes(cyc6)
#' ##-----------------------------------------------------------------
#' ##-----------------------------------------------------------------
#' ## 2 Repeats start. rep3,4
#' ## Backtrack to prev. level
#' 
enaCycle <- function (x) {

                                        #Initials

    if(class(x)!='network') {stop("x is not a network class object")}
    web <- as.matrix(x,attrname="flow")
    y <- x %v% "output"
    z <- x %v% "input"
    N <- length(y)

    TPTS <- apply(web,2,sum) + z

    F <- web/TPTS

    TST <- sum(web)+sum(y)+sum(z)
    df<-data.frame(NULL)
    df.cycle<-data.frame(0,0,'cycle', stringsAsFactors=FALSE)
#'##-----------------------------------------------------------------

                                        #Zero Global Variables
    NFST <- NEXNUM <- NCYC <- 0
    CYCS <- rep(0,N)

#'##-----------------------------------------------------------------

                                        #Start primary repeat loop
    repeat {
                                        # Zero all local variables
        NNEX  <- 0
        TCYCS <- rep(0,N)
        TMP   <- web*0
                                        #-----------------------------------

                                        #Count cycle arcs and determine exit from return
        NFWD <- NULL
        NTEMP <-NULL
        for (ii in 1:N) {#do 200 ii=1,N
            NFWD <- rep(0,N)
            NFWD[ii] <- 1
            for (k in 1:(N-1)) {
                for (i in 1:N) {
                    if(NFWD[i]>0) {next}
                    for (j in 1:N) {
                        if (NFWD[j] <1) {next}
                        if (web[j,i]<=0) {next}
                        NFWD[i] <- 1
                        break
                    }
                }
            }
            NTEMP[ii] <- 0
            for (i in 1:N) {
                if ((NFWD[i]>0) && (web[i,ii]>0)) {NTEMP[ii] <- NTEMP[ii]+1}
            }
        }#200
                                        # NTEMP will give the no. of Cycles ending in each node
        NSTP <- 0
        MAP <- NULL
        for (i in 1:N) {
            NMAX <- -1
            for (j in 1:N) {
                if (NTEMP[j]<=NMAX) {next}
                NMAX <- NTEMP[j]
                JMAX <- j
            }
            if (NMAX > 0) {NSTP <- NSTP + 1}
            NTEMP[JMAX] <- -2
            MAP[i] <- JMAX
        }
        #print(c('NSTP',NSTP))


        ### Condition for breaking/Exiting the primary repeat loop
        ###----------------------------------------
        if (NSTP<=0) {break} #breaks the primary repeat loop
        ###----------------------------------------

                                        # Start the NSTP2 While loop for Critical Arc determination
        ###----------------------------------------------------------------------------------------------
        NSTP2 <- 0
        repeat {
            slf.loop<-FALSE
            ARCMIN <- 10^25 #arbitrary min arc
            for (ir in 1:NSTP) {
                for (ic in 1:NSTP) {
                    IRTP <- MAP[ir]  # IRTP is the node to be searched at the ir'th position
                    ICTP <- MAP[ic]  # ICTP ----------------------------------ic'th --------
                    if (web[IRTP,ICTP] <= 0) {next}
                    if (web[IRTP,ICTP] >= ARCMIN) {next}
                    ARCMIN <- web[IRTP,ICTP]
                    IMIN <- IRTP
                    IM <- ir
                    JMIN <- ICTP
                    JM <- ic
                }
            }
            #print(ARCMIN)
            #print(min(web[web>0]))
            ### Exit from while(NSTP2<=0) if slf.loop
            if (IMIN == JMIN) {
                slf.loop <- TRUE
                break    #-----------------------BREAK THE WHILE repeat LOOP
            }
                                        #Make sure at least one cycle contains the current smallest arc web[IMIN,JMIN]
            NHALF <- (N/2)+1
            NFWD <- rep(0,N)
            NFWD[JMIN] <- 1
            ### find nodes from JMIN in fwd dirctn
            for (k in 1:NHALF) {
                for (i in 1:N) {
                    if (NFWD[i] >0) {next}
                    for (j in 1:N) {
                        if (NFWD[j] < 1) {next}
                        if (web[j,i] <= 0) {next}
                        NFWD[i] <- 1
                        break
                    }
                }
            }
            ### find nodes from IMIN in bckwd dirctn
            NODE <- rep(0,N)
            NODE[IMIN] <- 1
            for (k in 1:NHALF) {
                for (i in 1:N) {
                    if (NODE[i] > 0) {next}
                    for (j in 1:N) {
                        if (NODE[j] < 1) {next}
                        if (web[i,j] <= 0) {next}
                        NODE[i] <- 1
                        break
                    }
                }
            }
            ### find Common nodes aka members of the NEXUS
            NSTP2 <- 0  ###### ?????
            MAP2 <- rep(0,NSTP)
            for (i in 1:NSTP) {
                if ((NFWD[MAP[i]] <= 0) | (NODE[MAP[i]] <= 0)) {next}
                NSTP2 <- NSTP2 + 1
                MAP2[NSTP2] <- MAP[i]
            }
                                        # reorder mapping for IMIN and JMIN to come 1st and 2nd
            NFWD <- MAP2
            MAP2[1] <- IMIN
            MAP2[2] <- JMIN ##### SHORTER WAY POSSIBLE ####
            if (NSTP2 > 2) {
                INDX <- 2
                for (i in 1:NSTP2) {
                    if ((NFWD[i] == IMIN)||(NFWD[i]==JMIN)) {next}
                    INDX <- INDX+1
                    MAP2[INDX] <- NFWD[i]
                }
            }
            if (NSTP2 > 0) {break}
            web[IMIN,JMIN] = -web[IMIN,JMIN]
        } ###End of While Repeat NSTP2<=0
        ###----------------------------------------------------------------------------------------

                                        #IF the Critical arc is self loop
                                       #---------------------------------
        if(slf.loop == TRUE) {
            slf.loop <- FALSE
        	NCYC <- NCYC+1
            NNEX <- NNEX+1
            CYCS[1]<- CYCS[1]+web[IMIN,JMIN]
            WKARC<-F[IMIN,JMIN]*TPTS[IMIN]
            curr.slf.cyc<-noquote(c(NCYC,'.','(',IMIN,JMIN,')'))
            #print(curr.slf.cyc)
            #this.cycle <- rep(NA,N)
            this.cycle <- c(IMIN,JMIN)
            this.cycle <- paste(this.cycle, collapse='-')
            newcycle<-c(NCYC,(NEXNUM+1),this.cycle)
            df.cycle<-rbind(df.cycle,newcycle)

            web[IMIN,JMIN] <- 0
            NEXNUM <- NEXNUM+1
            curr.nexus <- c(NEXNUM,NNEX,IMIN,JMIN,WKARC)

            df<-rbind(df,curr.nexus) ### ('NEXUS', 'Cycles', weakarc, fromnode, tonode) ####################### df

            NFST <- 1
        }#End of if(slf.loop==TRUE)#

                                        #Begin Search for NEXUS defined by web[IMIN,JMIN] if not a self loop
        ###----------------------------------------------------------------------------------------

        else {
            WHOLE <- 0
                                        # Backtrack Routine Starts --------------
             # Backtrack Routine Starts --------------
            ### Initialize Node and Level
            LEVEL  <- 2
            NODE[1]<- 1
            NODE[2]<- 2
            skip.con.adv <- FALSE
            nex.com <- FALSE

            ### 2 Repeats start. rep1,2.
            repeat { #rep1
                repeat { #rep2
                    skip.con.chk <- FALSE

             ## Adv to nxt levels
                    if(skip.con.adv==FALSE) {
                        LM1         <- LEVEL
                        LEVEL       <- LEVEL+1
                        NODE[LEVEL] <- 1
                    }
#'## 2 Repeats start. rep3,4
                    repeat { #rep3
                        repeat { #rep4
                            ## Check for conn. b/w nodes at prsnt levels
                            if(skip.con.adv==FALSE){
                                NZ1  <- NODE[LM1]
                                KROW <- MAP2[NZ1]
                                NZ2  <- NODE[LEVEL]
                                KCOL <- MAP2[NZ2]
                                conn.chk <- FALSE
                                ## break rep4 and rep3 if conn exists
                                if(web[KROW,KCOL]>0 && skip.con.chk==FALSE) {
                                    conn.chk <- TRUE
                                    break #brk rep4
                                }
                            }
                            ## try next node in nxt level
                            NODE[LEVEL] <- NODE[LEVEL]+1
                            skip.con.adv <- FALSE
                            skip.con.chk <- FALSE
                            ## break rep4 after all levels are checked(NODE[level]>NSTP2)
                            if(NODE[LEVEL]>NSTP2) {break} #brk rep4
                        }#end of rep4
                        if(conn.chk==TRUE) {
                            conn.chk <- FALSE
                            break #rep3
                        }
#'## Backtrack to prev. level
                        LEVEL <- LEVEL-1
                        LM1   <- LEVEL-1
                        ## if further backtracking is impossible,
                                        #end search under weak arc. nexus complete
                        nex.com<-FALSE
                        if(LEVEL <= 2){
                            nex.com <- TRUE
                            break #break rep3
                        }
                        else {skip.con.chk<-TRUE} #goes to #420 to inc NODE[LEVEL]
                    }#end of rep3
                    if(nex.com==TRUE) {break} #break rep2
                    ##break rep2 if this conn completes cycle
                    if(NODE[LEVEL]==1) {break} #brk rep2
                    skip.con.adv <- FALSE
                    for (k in 1:LM1) {if (NODE[LEVEL] == NODE[k]) {skip.con.adv <- TRUE}}

                }#end of rep2
                if(nex.com==TRUE) {break} #brk rep1
                                        # -----------------BR Ends --------------
                                        # Calculate circuit prob
                WEIGHT <- 1
                for (kk in 1:LM1) {
                    KKP1 <- kk+1
                    KROW <- NODE[kk]
                    KCOL <- NODE[KKP1]
                    KROW <- MAP2[KROW]
                    KCOL <- MAP2[KCOL]
                    WEIGHT <- WEIGHT*F[KROW,KCOL]
                }
                                        # Add this weight
                for (kk in 1:LM1) {
                    KKP1 <- kk+1
                    KROW <- NODE[kk]
                    KCOL <- NODE[KKP1]
                    KROW <- MAP2[KROW]
                    KCOL <- MAP2[KCOL]
                    # Also, add this amount to the cycle distributions
                    TCYCS[LM1] <- TCYCS[LM1]+WEIGHT
                    TMP[KROW,KCOL] <- TMP[KROW,KCOL]+WEIGHT
                }
                                        # Report this cycle
                NNEX <- NNEX+1
                KTRY <- NNEX%%5000
                curr.prog <- noquote(c(NNEX,'Nexus cycles and Counting'))
                if(KTRY==0) {print(curr.prog)}
                NCYC <- NCYC+1
                if(NNEX>50){
                	skip.con.adv <- TRUE
                	next
                }
                L0 <- LM1+1
                for (kk in 1:L0) {
                    NTMP <- NODE[kk]
                    NTEMP[kk] <- MAP2[NTMP]
                }
                #curr.cycle <- noquote(c(NCYC,'.',NTEMP[1:L0]))
                #print(curr.cycle)
                this.cycle <- NTEMP[1:L0]
                #this.cycle[(L0+1):N]<-NA
                this.cycle <- paste(this.cycle, collapse='-')
                newcycle<-c(NCYC,(NEXNUM+1),this.cycle)
                df.cycle<-rbind(df.cycle,newcycle) #################------------------------------------df.cycle
                if(NNEX==50) {df.cycle<-rbind(df.cycle,rep('.',3))}
                skip.con.adv<-TRUE
            }#end of rep1
                                        #-----------------NEXUS COMPLETED---NEXUS REPEAT(rep1) ENDS HERE



                                        # Report this NEXUS
            WKARC=F[IMIN,JMIN]*TPTS[IMIN]
            NEXNUM=NEXNUM+1
            curr.nexus <- c(NEXNUM,NNEX,IMIN,JMIN,WKARC)
            df<-rbind(df,curr.nexus) ### ('NEXUS', 'Cycles', weakarc, fromnode, tonode) ####################### df
            #print(curr.nexus)
            # -------------------------------------
            # Normalize Probability Matrix & Subtract proper amounts from web
            PIVOT <- TMP[IMIN,JMIN]
            if(PIVOT<=0){print('Error in Normalizing Nexus Weights')}
            for(i in 1:N) {
                for(j in 1:N) {
                    if(web[i,j]<=0) {next}
                    web[i,j] <- web[i,j]-((TMP[i,j]/PIVOT)*ARCMIN)
                }
            }

            # Add proper amts to cycle distributions
            for(i in 1:N){CYCS[i]<-CYCS[i]+((TCYCS[i]/PIVOT)*ARCMIN)}               ############################## ERR.CHK CYCS depends on TCYCS, PIVOT, ARCMIN

            # Zero weak arc
            web[IMIN,JMIN] <- 0
            NFST <- 1
            if(WHOLE>1.00001) {print(c('Bad Sum Check = ',WHOLE))}


        }#End of else i.e. if not a slf.loop
    }#End of primary repeat loop

                                              # Report the Overall Results

    ### FIRST, UNCOVER ANY LINKS "HIDDEN" DURING SEARCH.
    web=abs(web)
    if(NFST!=0) {
        #cyc.rem <- noquote((c('A total of ',NCYC,'Cycles removed')))
        #print(cyc.rem)
        #print('Cycle Distributions')
        #print(CYCS)
        cycs<-CYCS
        CYC  <- sum(CYCS)
        CYCS <- (CYCS/TST)
        #print('Normalized Distribution')
        #print(CYCS)
        TEMP <- CYC/TST
        #print(c('cycling index is',TEMP))
        ResidualFlows<-web
        AggregatedCycles<-(as.matrix(x, attrname = 'flow')) - ResidualFlows
        colnames(df)<-c('NEXUS', 'CYCLES','W.arc.From','W.arc.To', 'W.arc.Flow')
        #colnames(df.cycle)<-rep(' ',(N+2))
        colnames(df.cycle)<-c('CYCLE','NEXUS','NODES')
        #df.cycle[is.na(df.cycle)==TRUE]<- ' '
        df.cycle<-df.cycle[-1,]
        rw<-row.names(df.cycle); rw<-as.numeric(rw)
        row.names(df.cycle) <- 2:max(rw)-1
        NCYCS<-NCYC; NNEX<-NEXNUM; CI<-TEMP
        ns <- cbind(NCYCS, NNEX, CI)
        out <- list(Table.cycle=df.cycle,Table.nexus=df,CycleDist = cycs, NormDist=CYCS, ResidualFlows=web, AggregatedCycles=AggregatedCycles, ns=ns)
        return(out)
    }#end of if (NFST!=0)
    else {
        NCYCS<-NCYC;NNEX<-NEXNUM; CI<-0
        ns <- cbind(NCYCS, NNEX, CI)
        out <- list(ResidualFlowks=web,ns=ns)
        return(out)
      }


}#END OF FUNCTION
#' environ --- conducts environ analysis 
#' INPUT = network object
#' OUTPUT = input and/or output environs
#' 
#' M. Lau July 2011 | DEH edited Feb 2013
#' ---------------------------------------------------







#' environ --- conducts environ analysis INPUT = network object OUTPUT = input
#' and/or output environs
#' 
#' M. Lau July 2011 | DEH edited Feb 2013
#' --------------------------------------------------- environ --- conducts
#' environ analysis INPUT = network object OUTPUT = input and/or output
#' environs
#' 
#' M. Lau July 2011 | DEH edited Feb 2013
#' --------------------------------------------------- Ecological Network
#' Environs
#' 
#' Calculates the environs for an ecological network.
#' 
#' @param x A network object.
#' @param input Should the input environ be calculated?
#' @param output Should the output environ be calculated?
#' @param type Specifies the type of environs ("unit" or "realized") to be
#' calculated.
#' @param err.tol Error threshold for numerical error fluctuations in flows.
#' Values below err.tol will be set to zero.
#' @param balance.override Logical specifying whether (TRUE) or not (FALSE) the
#' model needs to be balanced prior to calculations. If TRUE and the model is
#' not balanced, environs will not be calculated.
#' @return The function returns the input, output or both environs depending
#' upon which were requested.
#' @author Stuart R. Borrett Matthew K. Lau
#' @references Fath, B.D. and S.R. Borrett. 2006. A MATLAB function for network
#' environ analysis. Environmental Modelling & Software 21:375-405.
#' @examples
#' 
#' 
#' 
#' data(troModels)
#' enaEnviron(troModels[[6]])
#' 
#' 
#' 
#' @export enaEnviron
enaEnviron <- function(x,input=TRUE,output=TRUE,type='unit',err.tol=1e-10,balance.override=FALSE){
                                        #check for network class
  if (class(x) != 'network'){warning('x is not a network class object')}
                                        #Check for balancing
  if (balance.override){}else{
    if (any(list.network.attributes(x) == 'balanced') == FALSE){x%n%'balanced' <- ssCheck(x)}
    if (x%n%'balanced' == FALSE){warning('Model is not balanced'); stop}
  }
                                        #Assume 'rc' orientation of flows
                                        #Don't transpose flows for calculations, they will be transposed in enaFlow
                                        #calculate enaFlow with RC input
  user.orient <- get.orient()
  set.orient('rc')
  Flow <- enaFlow(x)
  set.orient(user.orient)
                                        #Transpose for calculations in Patten school
  Flow$G <- t(Flow$G)
  Flow$GP <- t(Flow$GP)
  Flow$N <- t(Flow$N)
  Flow$NP <- t(Flow$NP)
                                        #Unit environ calculations
  if (input){
                                        #Input perspective
    EP <- list()
    for (i in (1:nrow(Flow$NP))){
      dNP <- diag(Flow$NP[i,]) #diagonalized N matrix
      EP[[i]] <- dNP %*% Flow$GP         # calculate internal environ flows
      EP[[i]] <- EP[[i]] - dNP      # place negative environ throughflows on the principle diagonal
      EP[[i]] <- cbind(EP[[i]],apply((-EP[[i]]),1,sum)) #attach z column
      EP[[i]] <- rbind(EP[[i]],c(apply((-EP[[i]]),2,sum))) #attach y row
      EP[[i]][nrow(EP[[i]]),ncol(EP[[i]])] <- 0 #add zero to bottom right corner to complete matrix
      EP[[i]][abs(EP[[i]]) < err.tol] <- 0 #ignore numerical error
                                        #add labels to matrices
      labels <- c(rownames(Flow$GP))
      colnames(EP[[i]]) = c(labels, 'z')
      rownames(EP[[i]]) = c(labels, 'y')
    }
                                        #add environ names
    names(EP) <- labels
  }
  if (output){
                                        #Output perspective
    E <- list()
    for (i in (1:nrow(Flow$N))){
      dN <- diag(Flow$N[,i]) #diagonalized N matrix
      E[[i]] <- Flow$G %*% dN
      E[[i]] <- E[[i]] - dN
      E[[i]] <- cbind(E[[i]], apply((-E[[i]]), 1, sum))
      E[[i]] <- rbind(E[[i]], c(apply((-E[[i]]), 2, sum)))
      E[[i]][nrow(E[[i]]), ncol(E[[i]])] <- 0 
      E[[i]][abs(E[[i]]) < err.tol] <- 0 
                                        #add labels to matrices
      labels <- c(rownames(Flow$G))
      colnames(E[[i]]) = c(labels, 'z')
      rownames(E[[i]]) = c(labels, 'y')
    }
                                        #add environ names
    names(E) <- labels
  }
                                        #Realized environ calculations
  if (type == 'realized'){
                                        #Input perspective
    if (input){
      for (i in (1:nrow(Flow$N))){
        EP[[i]] <- EP[[i]]*unpack(x)$y[i] #Construct realized environ
      }
    }
                                        #Output perspective
    if (output){
      for (i in (1:nrow(Flow$N))){
        E[[i]] <- E[[i]]*unpack(x)$z[i]
      }
    }
  }
                                        #Wrap-up output into list  
  if (input & output){
    out <- list('input' = EP,'output' = E)
  } else if (input & output == FALSE){
    out <- EP
  } else if (input == FALSE & output){
    out <- E
  }      
  
  if (type != 'unit' && type!= 'realized'){
    print('WARNING: Invalid input in type, input ignored')
  }
                                        #re-orient matrices
  if (user.orient == 'rc'){
    for (i in 1:length(out)){
      for (j in 1:length(out[[i]])){
        out[[i]][[j]] <- t(out[[i]][[j]])
      }
    }
  }else{}
                                        #output
  return(out)
}

#' enaFlow --- flow analysis
#' INPUT = network object
#' OUTPUT = list of flow statistics
#'
#' M. Lau | July 2011
#' ---------------------------------------------------







#' enaFlow --- flow analysis INPUT = network object OUTPUT = list of flow
#' statistics
#' 
#' M. Lau | July 2011 ---------------------------------------------------
#' enaFlow --- flow analysis INPUT = network object OUTPUT = list of flow
#' statistics
#' 
#' M. Lau | July 2011 --------------------------------------------------- Flow
#' Analyses of Ecological Networks
#' 
#' Performs the primary throughflow analysis developed for input-output
#' systems.  It returns a vector of throughflows, the input and output oriented
#' matrices for "direct flow intensities" and "integral flow intensities", and
#' a set of flow based network statistics.  Included in the network statistics
#' are a set of measures that describe the diversity of flows in the ecosystem
#' (Ulanowicz's Ascendendy measures).
#' 
#' @param x a network object.  This includes all weighted flows into and out of
#' each node.
#' @param zero.na LOGICAL: should NA values be converted to zeros.
#' @param balance.override Flow analysis assumes the network model is at
#' steady-state (inputs = outputs).  Setting balance.override = TRUE allows the
#' function to be run on unbalanced models.
#' @return \item{T}{vector of node throughflows - total amount of energy-matter
#' flowing into or out of each node} \item{G}{matrix of the output oriented
#' direct flow intensities} \item{GP}{matrix of the input oriented direct flow
#' intensities} \item{N}{matrix of the ouput oriented integral
#' (boundary+direct+indirect) flow intensities} \item{NP}{matrix of the input
#' oriented integral flow intensities} \item{ns}{vector of flow based network
#' statistics.  These include "Boundary" the total input into or output from
#' the system, "TST" the total system throughflow, "TSTp" total system
#' throughPUT,"APL" is the network aggradation TST/Boundary which is also
#' called average path length, "FCI" (Finn Cycling Index) is a metric of the
#' amount of cycling in a system, "BFI" is the boundary flow intensity
#' Boundary/TST, "DFI" is the direct flow intensity Direct/TST, "IFI" is the
#' indirect flow intensity Indirect/TST, "ID.F" is the realized indirect to
#' direct flow intensity, "ID.F.I" is the input idealized indirect flow
#' intensity, "id.F.O"is the output idealized indirect flow intensity, "HMG.I"
#' is the input network homogenization, "HMG.O" is the output network
#' homogenization, "AMP.I" is the strong measure of input network amplifiation,
#' "AMP.O" is the strong measure of output network amplification, "mode0.F" is
#' the boundary flow - flow that reaches a compartment from across the system
#' boundary, "mode1.F" is internal first passage flow, "mode2.F" is cycled
#' flow, "mode3.F" is the dissipative eqivalent to mode2, and "mode4.F" is the
#' dissipative equivalent ot mode0. "H" is the total flow diversity (MacArthur
#' 1955).  Uses the Shannon Information measure (aka Boltzmann entropy) applied
#' to the individual flows, "AMI", is the Average Mutual Information (AMI) in a
#' network. "Hr", is the residual uncertainty that remains about the flow
#' distribution once the ecosystem structure is specified (Hr = H - AMI),
#' "ASC", Returns the ascendnecy of a network, "OH", is the overhead of a
#' network (Ulanowicz 2004), "CAP", is the capacity of a network, "ACS.CAP", is
#' the proportion of capacity used by ascendency, "OH.CAP", Returns the
#' proportion of capacity used by overhead, "robustness", is the robustness of
#' the network, "ELD" Returns the Effective Link Density of the network(c)
#' (Ulanowicz et al. 2014), "TD", Returns the Trophic Depth of the network(r)
#' (Ulanowicz et al. 2014), "A.input", Returns the input ascendnecy of a
#' network, "A.internal", Returns the internal ascendnecy of a network,
#' "A.export", Returns the export ascendnecy of a network, "A.respiration",
#' Returns the respiration ascendnecy of a network, "OH.input", Returns the
#' input overhead of a network, "OH.internal", Returns the internal overhead of
#' a network, "OH.export", Returns the export overhead of a network,
#' "OH.respiration", Returns the respiration overhead of a network,
#' "CAP.input", Returns the input capacity of a network, "CAP.internal",
#' Returns the internal capacity of a network, "CAP.export", Returns the export
#' capacity of a network, "CAP.respiration", Returns the respiration capacity
#' of a network}
#' @author Matthew K. Lau Stuart R. Borrett
#' @seealso
#' \code{\link{read.scor},\link{read.wand},\link{enaStorage},\link{enaUtility}}
#' @references Borrett, S. R., Freeze, M. A., 2011. Reconnecting environs to
#' their environment. Ecol. Model. 222, 2393-2403.
#' 
#' Fath, B. D., Borrett, S. R. 2006. A Matlab function for Network Environ
#' Analysis.  Environ. Model. Softw. 21, 375-405.
#' 
#' Fath, B. D., Patten, B. C., 1999. Review of the foundations of network
#' environ analysis. Ecosystems 2, 167-179.
#' 
#' Finn, J. T., 1976. Measures of ecosystem structure and function derived from
#' analysis of flows. J. Theor. Biol. 56, 363-380.
#' 
#' Patten, B.C. Higashi, M., Burns, T. P. 1990. Trophic dynamics in ecosystem
#' networks: significance of cycles and storage.  Ecol. Model. 51, 1-28.
#' 
#' Schramski, J. R., Kazanci, C., Tollner, E. W., 2011. Network environ theory,
#' simulation and EcoNet 2.0. Environ. Model. Softw. 26, 419-428.
#' 
#' Ulanowicz, R.E., 2004. Quantitative methods for ecological network analysis.
#' Comput. Biol. Chem. 28, 321-33
#' 
#' Ulanowicz, R.E., Holt, R.D., Barfield, M., 2014. Limits on ecosystem trophic
#' complexity: insights from ecological network analysis.  Ecology Letters
#' 17:127-136.
#' @examples
#' 
#' 
#' 
#' data(troModels)
#' F = enaFlow(troModels[[6]])  # completes the full analysis
#' F$ns  # returns just the network statisics
#' 
#' 
#' 
#' @export enaFlow
enaFlow <- function(x,zero.na=TRUE,balance.override=FALSE){
                                        #Check for network class
  if (class(x) != 'network'){warning('x is not a network class object')}

                                        #Check for balancing
  if (balance.override){}else{
    if (any(list.network.attributes(x) == 'balanced') == FALSE){x%n%'balanced' <- ssCheck(x)}
    if (x%n%'balanced' == FALSE){warning('Model is not balanced'); stop}
  }

                                        # unpack model
  Flow <- t(as.matrix(x,attrname = 'flow')) #flows
  input <- x%v%'input' #inputs
  stor <- x%v%'storage' #storage values

  n <- nrow(Flow)      # number of nodes
  I <- diag(1,nrow(Flow),ncol(Flow))          # create identity matrix
  T. <- apply(Flow,1,sum) + input;   # input throughflow (assuming steady state)

                                        #compute the intercompartmental flows
  GP <- Flow / T.  #Input perspective
  G <- t(t(Flow) / T.)  #Output perspective

                                        #check and replace NA values with 0 if zero.na
  if (zero.na){
    GP[is.na(GP)] <- 0
    G[is.na(G)] <- 0
    GP[is.infinite(GP)] <- 0
    G[is.infinite(G)] <- 0
  }

                                        #compute the integral flows
  NP <- ginv((I - GP))
  rownames(NP) <- colnames(NP) <- colnames(GP)
  N <- ginv((I - G))
  rownames(N) <- colnames(N) <- colnames(G)

  ## the ginv function creates noticible numeric error.  I am removing some of it here by rounding
  tol <- 10
  N <- round(N,tol)
  NP <- round(NP,tol)


  ## Network Statistics
  TST <- sum(T.)  # total system throughflow
  TSTp <- sum(Flow) + sum(x%v%'input') + sum(x%v%'output') # total system throughput

  Boundary <- sum(input)
  APL <- TST/Boundary  # Average Path Lenght (Finn 1976; aka network
                      # aggradation, multiplier effect)

                                        # Finn Cycling Index
  p <- as.matrix(rep(1,n),nrow=n)
  dN <- diag(N)
  TSTc <- sum((dN-p)/dN *T.)
  FCI <- TSTc/TST

                                        # non-locality (realized)
  direct <- sum(G %*% input)
  indirect <- sum((N - I - G) %*% input)
  ID.F <- indirect/direct
  BFI <- Boundary/TST
  DFI <- sum(G %*% input) / TST
  IFI <- indirect/TST
                                        # non-locality (idealized)
  ID.F.O <- sum(N-I-G)/sum(G)
  ID.F.I <- sum(NP-I-GP)/sum(GP)
                                        # HMG
  HMG.O <- ( sd(as.vector(G)) / mean(G) ) / ( sd(as.vector(N)) / mean(N) )
  HMG.I <- ( sd(as.vector(GP)) / mean(GP) ) / ( sd(as.vector(NP)) / mean(NP) )
                                        # Amplification
  AMP.O <- length(which( (N - diag(diag(N))) > 1))
  AMP.I <- length(which( (NP - diag(diag(NP))) > 1))
                                        # MODE ANALYSIS
                                        # This is built from Fath's original MODE program
  mode0.F <- Boundary                     # boundary flow
  mode1.F <- sum(ginv(diag(diag(N))) %*% N %*% diag(input) - diag(as.vector(I%*%input))) # internal first passage flow
  mode2.F <- sum((diag(diag(N))-I) %*% ginv(diag(diag(N))) %*% N %*% diag(input))  # cycled flow
  mode3.F <- mode1.F                      # dissipative equivalent to mode 1
  mode4.F <- mode0.F                    # dissipative equivalent to mode 0
                                        #re-orientation
  orient <- get.orient()
  if (orient == 'rc'){
      G <- t(G)
      GP <- t(GP)
      N <- t(N)
      NP <- t(NP)
  }


  asc <- enaAscendency(x)
                                        #network statistics
  ns <- cbind(Boundary,TST,TSTp,APL,FCI,
              BFI,DFI,IFI,
              ID.F,ID.F.I,ID.F.O,
              HMG.I,HMG.O,
              AMP.I,AMP.O,
              mode0.F,mode1.F,mode2.F,mode3.F,mode4.F, asc)

                                        #output
  return(list('T'=T.,'G'=G,'GP'=GP,'N'=N,'NP'=NP,'ns'=ns))
}


#' enaMTI --- Mixed Trophic Impacts Analysis
#' follows Ulanowicz and Puccia, 1990.
#' INPUT = network object
#' OUTPUT = list of trophic impact statistics
#' Borrett | June 2012, MKL | July 2013
#' ------------------------------------







#' enaMTI --- Mixed Trophic Impacts Analysis follows Ulanowicz and Puccia,
#' 1990. INPUT = network object OUTPUT = list of trophic impact statistics
#' Borrett | June 2012, MKL | July 2013 ------------------------------------
#' enaMTI --- Mixed Trophic Impacts Analysis follows Ulanowicz and Puccia,
#' 1990. INPUT = network object OUTPUT = list of trophic impact statistics
#' Borrett | June 2012, MKL | July 2013 ------------------------------------
#' Mixed Trophic Impacts (MTI) Analysis
#' 
#' Calculates the Mixed Trophic Impacts of one species on another in the given
#' ecosystem model following the algorithm of Ulanowicz and Puccia (1990). This
#' considers both the direct and indirect trophic impacts.
#' 
#' 
#' @param x a network object.  This includes all weighte dflows into and out of
#' each node.  It must also include the "Living" vector that identifies the
#' living (TRUE/FALSE) status of each node.
#' @param eigen.check LOGICAL: should the dominant eigen value be checked?  By
#' default, the function will not return utility values if the eigenvalue is
#' larger than one; however, if eigen.check is set to FALSE, then the function
#' will be applied regardless of the mathematic concern.
#' @param zero.na A logical parameter that specifies if NAs generated in the
#' analysis should be reset to zero.  The default is TRUE.
#' @param balance.override Mixed Trophic Impacts analysis builds on flow
#' analysis and thus assumes the network model is at steady-state (inputs =
#' outputs).  Setting balance.override = TRUE allows the function to be run on
#' unbalanced models, though this is unadvised.
#' @return \item{G}{output-oriented direct flow intensity matrix as in enaFlow,
#' except oriented from row to column.} \item{FP}{input-oriented direct flow
#' intensity matrix similar to enaFlow; however, the calculation exclude
#' respiration losses from the throughflow in the denominator to focus on NET
#' production.  Also, if the receiver compartment is not living, the flux
#' intensity is set to zero.} \item{Q}{direct net trophic impacts (G-t(FP)).}
#' \item{M}{Total (direct and indirect) tropic impacts of compartment i on j.}
#' @note This and other Ulanowicz school functions require that export and
#' respiration components of output be separately quantified.
#' 
#' This analysis is similar in concept to the ENA Utility analysis.
#' 
#' With regard to the eigen.check argument, like enaFlow, enaStorage and
#' enaUtility, this analysis considers the trophic impact propigated over path
#' lengths ranging for zero to infinity.  For the analysis to work properly,
#' the path sequence must converge.  This function checks to see if the path
#' sequence is convergent by finding the dominant eigenvalue of the direct
#' matrix.  If this eigenvalue is less than 1, the sequence is convergent and
#' the analysis can be applied; if the dominant eigenvalue is greater than one,
#' then the anlysis cannot be applied.
#' @author Stuart R. Borrett Matthew K. Lau
#' @seealso \code{\link{enaFlow},\link{enaUtility}}
#' @references %% ~put references to the literature/web site here ~ Ulanowicz,
#' R.E. and C.J. Puccia.  1990. Mixed trophic impacts in ecosystems.  Coenoses
#' 5, 7--16.
#' @examples
#' 
#' 
#' 
#' data(troModels)
#' mti <- enaMTI(troModels[[6]])
#' attributes(mti)
#' 
#' 
#' 
#' @export enaMTI
enaMTI <- function(x,eigen.check=TRUE,zero.na=TRUE, balance.override=FALSE){
                                        #Check for network class
  if (class(x) != 'network'){warning('x is not a network class object')}
                                        #Data checks
  if (any(is.na(x%v%'respiration'))){
    G <- FP <- Q <- M <- as.matrix(x, attrname = 'flow')
    G[is.na(G)==FALSE] <- FP[is.na(FP)==FALSE] <- Q[is.na(Q)==FALSE] <- M[is.na(M)==FALSE] <- NA
    out <- list('G'=G,'FP'=FP,'Q'=Q,'M'=M)
    warning('Model is missing respiration. Output is NA.')
  }else{
                                        #Check for balancing
    if (balance.override){}else{
      if (any(list.network.attributes(x) == 'balanced') == FALSE){x%n%'balanced' <- ssCheck(x)}
      if (x%n%'balanced' == FALSE){warning('Model is not balanced'); stop}
    }
                                        #Unpack
    Flow <- as.matrix(x, attrname = 'flow')  #flows
    input <- x%v%'input' #inputs
    output <- x%v%'output'
    resp <- x%v%'respiration'
    stor <- x%v%'storage' #storage values
    I <- Flow*0;
    diag(I) <- 1 #create the identity matrix
    T. <- input + apply(Flow,2,sum)
                                        #
    G <- t(t(Flow)/T.)        # input oriented direct flow intensity matrix
    FP <- Flow / (T.-resp)    # modified output oriented direct flow intensity matrix.  Authors exclude respiration to divide only by the NET production of the compartment.

                                        # check and replace NA values with 0 if zero.na
    if (zero.na){
      G[is.na(G)] <- 0
      FP[is.na(FP)] <- 0
    }

    # Make infinity values equal to zero
    G[is.infinite(G)] <- 0
    FP[is.infinite(FP)] <- 0

    # Set FP to zero when receiver compartment (j) is non-living
    FP[,which(x%v%'living'==FALSE)] <- 0
    Q <- G - t(FP)
    dom1Q <- abs(eigen(Q)$values[1])
    if(dom1Q <= 1 ){
      M <- ginv(I-Q)-I              # Total Impacts of i on j.
    } else {
      if(eigen.check==FALSE){
        M <- ginv(I-Q)-I              # Total Impacts of i on j.
      } else { M <- NA}
  }

    if(!any(is.na(M))){
        r <- relationalChange(Q,M)
        IR <- r$Integral.Relations
        r.table <- r$Relations.Table
        names(r.table) <- c("From","To","Net (direct)","Mixed (integral)","changed")
        rownames(r.table) <- c(1:dim(r.table)[1])
    } else {
        IR <- NA
        r.table <- NA
    }

    out <- list('G'=G,'FP'=FP,'Q'=Q,'M'=M,
                "Integral.Relations" = IR,
                "Relational.Table"=r.table)
  }
    return(out)
}


#' Bigeochemical Cycling Models
#' 
#' A set of 43 biogeochemical cycling models compiled by the SEE Lab at UNCW.
#' 
#' 
#' @name bgcModels
#' @docType data
#' @references Borrett, S. R., and M. K. Lau. In Prep. enaR: An R package for
#' Ecological Network Analysis. Ecological Modeling and Software.
#' @keywords datasets
NULL





#' Ecosystem Model Information
#' 
#' Model information for the set of ecosystem models compiled by the SEE Lab at
#' UNCW.
#' 
#' 
#' @name enaModelInfo
#' @docType data
#' @references Borrett, S. R., and M. K. Lau. In Prep. enaR: An R package for
#' Ecological Network Analysis. Ecological Modeling and Software.
#' @keywords datasets
NULL





#' Ecosystem Models
#' 
#' A set of ecosystem models compiled by the SEE Lab at UNCW.
#' 
#' 
#' @name enaModels
#' @docType data
#' @references Borrett, S. R., and M. K. Lau. In Prep. enaR: An R package for
#' Ecological Network Analysis. Ecological Modeling and Software.
#' @keywords datasets
NULL





#' Tools for Ecological Network Analysis (ena)
#' 
#' This package compiles functions for the analysis of ecological networks,
#' building on tools previously developed in the MatLab language Borrett 2006)
#' with multiple additions of functionality.
#' 
#' \tabular{ll}{ Package: \tab enaR \cr Type: \tab Package\cr Version: \tab
#' 2.9.3\cr Date: \tab 2015-12-09\cr License: \tab GPL-3\cr }
#' 
#' @name enaR-package
#' @aliases enaR enaR-package
#' @docType package
#' @author Authors: Stuart R. Borrett, Matthew K. Lau, Pawandeep Singh, David
#' E. Hines Maintainer: Matthew K. Lau <enaR.maintainer@@gmail.com>
#' @seealso % ~~ Optional links to other man pages, e.g. ~~ % ~~
#' \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~ %
#' \code{\link[network:statnet-package]{statnet}}
#' %\code{\link[igraph:igraph-package]{igraph}}
#' \code{\link[network:network-package]{network}}
#' @references Borrett SR and Lau MK 2014. enaR: An r package for Ecosystem
#' Network Analysis. Methods in Ecology and Evolution 5:1206-1213.
#' @keywords package
NULL





#' Sub-set of the Larger Ecosystem Models
#' 
#' A set of ecosystem models compiled by the SEE Lab at UNCW.
#' 
#' 
#' @name m.list
#' @docType data
#' @references Borrett, S. R., and M. K. Lau. In Prep. enaR: An R package for
#' Ecological Network Analysis. Ecological Modeling and Software.
#' @keywords datasets
NULL





#' Intertidal Oyster Reef Ecosystem Model
#' 
#' Intertidal oyster reef ecosystem model created by Dame and Patten (1981).
#' Data were taken from Patten (1985).  Model flows are in kcal m^-2 day^-1;
#' storage data is kcal m^-2.
#' 
#' 
#' @name oyster
#' @docType data
#' @references Dame, R. F., and B. C. Patten. 1981. Analysis of energy flows in
#' an intertidal oyster reef. Marine Ecology Progress Series 5:115-124.
#' 
#' Patten, B. C. 1985. Energy cycling, length of food chains, and direct versus
#' indirect effects in ecosystems. Can. Bull. Fish. Aqu. Sci. 213:119-138.
#' @keywords datasets
NULL





#' Trophic Models
#' 
#' A set of 58 trophic models compiled by the SEE Lab at UNCW.
#' 
#' 
#' @name troModels
#' @docType data
#' @references Borrett, S. R., and M. K. Lau. In Prep. enaR: An R package for
#' Ecological Network Analysis. Ecological Modeling and Software.
#' @keywords datasets
NULL



#' enaStorage --- storage analysis
#' INPUT = network object
#' OUTPUT = list of storage statistics
#'
#' M. Lau | July 2011
#' ---------------------------------------------------







#' enaStorage --- storage analysis INPUT = network object OUTPUT = list of
#' storage statistics
#' 
#' M. Lau | July 2011 ---------------------------------------------------
#' enaStorage --- storage analysis INPUT = network object OUTPUT = list of
#' storage statistics
#' 
#' M. Lau | July 2011 ---------------------------------------------------
#' Storage Analyses of Ecological Networks
#' 
#' Calculates storage-based Ecological Network Analyses.
#' 
#' @param x A network object.  This This includes all weighted flows into and
#' out of each vertex as well as the amount of energy--matter stored at each
#' vertex.
#' @param balance.override LOGICAL: should an imbalanced model be analyzed?  If
#' FALSE, the functions checks to make sure the network model provided is at
#' steady-state.  If TRUE, then the function will run without ensuring that the
#' model meets the steady-state assumption.
#' @return \item{X}{The storage values themselves.} \item{C}{output or
#' donor-storage normalized output-oriented direct flow intensity matrix
#' (Jacobian community matrix)} \item{S}{dimensionalized integral output
#' community matrix} \item{Q}{integral output storage matrix - non-dimensional}
#' \item{CP}{input or recipient-storage normalized oriented flow intensity
#' matrix (Jacobian community matrix)} \item{SP}{dimensionalized integral input
#' community matrix} \item{QP}{integral input storage matrix - non-dimensional}
#' \item{dt}{selected time step to create P, PP, Q and QP - smallest whole
#' number to make diag(C) nonnegative} \item{ns}{vector of the storage based
#' whole system network statistics.  These statistics include total system
#' storage (TSS), storage cycling index (CIS), Boundary storage intensity
#' (BSI), Direct storage intensity (DSI), Indirect storage intensity (ISI),
#' realized ratio of indirect-to-direct storage (ID.S), unit input-oriented
#' ratio of indirect-to-direct storage intensities (IDS.I), unit output ratio
#' of indirect-to-direct storage intensities (IDS.O), input-oriented
#' storage-based network homogenization (HMG.S.I), output-oriented
#' storage-based network homogenization (HMG.S.O), input-oriented storage-based
#' network amplification (AMP.S.I), output-oriented storage-based network
#' amplification (AMP.S.O), Storage from Boundary flow (mode0.S), storage from
#' internal first passage flow (mode1.S), storage from cycled flow (mode2.S),
#' dissipative equivalent to mode1.S (mode3.S), dissipative equivalent to
#' mode0.S (mode4.S).}
#' @author Matthew K. Lau Stuart R. Borrett
#' @seealso
#' \code{\link{read.scor},\link{read.wand},\link{enaFlow},\link{enaUtility}}
#' @references Matis, J. H., Patten, B. C. 1981. Environ analysis of linear
#' compartmental systems: the static, time invariant case.  Bulletin of the
#' International Statistical Institute, 48: 527-565.
#' 
#' Fath, B. D., Patten, B. C. 1999.  Review of the foundations of network
#' enviorn analysis.  Ecosystems 2:167-179.
#' 
#' Fath, B. D. Patten, B. C., Choi, J. 2001.  Compementarity of ecological goal
#' functions.  Journal of Theoretical Biology 208: 493-506.
#' 
#' Fath, B. D., Borrett, S. R. 2006. A MATLAB function for Network Environ
#' Analysis.  Environmental Modelling & Software 21:375-405
#' @keywords enaFlow read.scor
#' @export enaStorage
#' @examples
#' data(oyster)
#' S <- enaStorage(oyster)
#' attributes(S)
enaStorage <- function(x,balance.override=FALSE){
                                        #Missing Data Check
  if (any(is.na(x%v%'storage'))){
    warning('This function requires quantified storage values.')
  }else{
                                        #Check for network class
    if (class(x) != 'network'){warning('x is not a network class object')}

                                        #Check for balancing
    if (balance.override){}else{
      if (any(list.network.attributes(x) == 'balanced') == FALSE){x%n%'balanced' <- ssCheck(x)}
      if (x%n%'balanced'){}else{stop('Model is not balanced')}
    }
                                        #unpack data from x
    Flow <- t(as.matrix(x, attrname = 'flow'))  #flows
                                        #continue unpacking
    input <- x%v%'input' #inputs
    stor <- x%v%'storage' #storage values
    T. <- apply(Flow,1,sum) + input
    FD <- Flow - diag(T.) #flow matrix with negative throughflows on the diagonal
    I <- diag(1,nrow(Flow),ncol(Flow)) #create the identity matrix

                                        #Compute the Jacobian matrix
    C <- FD %*% ginv(diag(stor)) #output matrix
    CP <- ginv(diag(stor)) %*% FD #input matrix

                                        #smallest whole number to make diag(C) nonnegative
    dt <- -1 / floor(min(diag(C)))

                                        #calculating the storage-specific, output-oriented, intercompartmental flows (P)
    P <- I + C*dt
    PP <- I + CP*dt
                                        #calculating the dimensionalized integral output and input matrices -- expected residence times (Barber 1979)
    S <- -ginv(C) #output
    SP <- -ginv(CP) #input

    tol <- 10
    S <- round(S,10)
    SP <- round(SP,10)

    #' calculate variance of expected residence times (Barber 1979)
    VS <- 2 * ( t(S) %*% diag(diag(S))) - t(S)^2
    VSP <- 2 * (t(SP) %*% diag(diag(SP))) - t(SP)^2


                                        #calculating the integral storage intensity matrix (Q)
    Q <- ginv(I - P) #output
    QP <- ginv(I - PP) #input

  ## the ginv function creates noticible numeric error.  I am removing some of it here by rounding
    tol <- 10
    Q <- round(Q,tol)
    QP <- round(QP,tol)


    dQ <- diag(Q) #diagonal of integral output storage matrix which is the same for input (i.e. diag(QP))

                                        #naming row and columns
    rownames(C) <- colnames(C) <- rownames(Flow)
    rownames(CP) <- colnames(CP) <- rownames(Flow)
    rownames(P) <- colnames(P) <- rownames(Flow)
    rownames(S) <- colnames(S) <- rownames(Flow)
    rownames(VS) <- colnames(VS) <- rownames(Flow)
    rownames(Q) <- colnames(Q) <- rownames(Flow)
    rownames(CP) <- colnames(CP) <- rownames(Flow)
    rownames(PP) <- colnames(PP) <- rownames(Flow)
    rownames(SP) <- colnames(SP) <- rownames(Flow)
    rownames(VSP) <- colnames(VSP) <- rownames(Flow)
    rownames(QP) <- colnames(QP) <- rownames(Flow)

    ##Storage Environ Properties
                                        #eigen analysis
    e <- eigen(P)$values
    lam1P <- e[1]
    rhoP <- e[1] / e[2]

    eP <- eigen(PP)$values
    lam1PP <- eP[1]
    rhoPP <- eP[1] / eP[2]

    TSS <- sum(stor) #total system storage
    TSScs <- sum(((dQ-1)/dQ)%*%stor) #cycled (mode 2) storage
    CIS <- TSScs / TSS #cucling index (storage)

    #' Amplification parameter
    NAS <- length((Q-diag(diag(Q)))[(Q-diag(diag(Q))) > 1])
    NASP <- length((QP-diag(diag(QP)))[(QP-diag(diag(QP))) > 1])

    #' Indirect effects parameter (srb fix 8.3.2011)
    ID.S.O <- sum(Q-I-P) /sum(P)
    ID.S.I <- sum(QP-I-PP) / sum(PP) #indirect to direct ratio (input matrix)

    #' Indirect effects parameter (realized)  (srb fix 8.3.2011)
    ID.S <- sum(dt*(as.matrix((Q-I-P)) %*% input)) / sum(dt*as.matrix(P)%*% input ) #indirect to direct ratio (output matrix)

    #' Tripartite walk-length division of storage
    BSI = sum( I %*% input *dt) / TSS
    DSI = sum( P %*% input *dt) / TSS
    ISI = sum( (Q-I-P) %*% input *dt) / TSS

    #' Homogenization parameter
    CVP <- sd(as.numeric(P)) / mean(P) #Coefficient of variation for G
    CVQ <- sd(as.numeric(Q)) / mean(Q)  #Coefficient of variation for N
    HMG.S.O <- CVP / CVQ #homogenization parameter (output storage)

    CVPP <- sd(as.numeric(PP)) / mean(PP) #Coefficient of variation for GP
    CVQP <- sd(as.numeric(QP)) / mean(QP) #Coefficient of variation for NP
    HMG.S.I <- CVPP / CVQP #homogenization paraemeter (input storage)

    #' Network Aggradation
    AGG.S <- TSS / sum(input) #network aggradation -- average amount of storage per system input

    #' MODE Partition (Fath et al. 2001)
    z <- unpack(x)$z
    mode0.S = sum(z*dt)  # storage from boundary input flow
    mode1.S = sum( (ginv(diag(diag(Q))) %*% Q - I) %*% diag(z) * dt ) # storage from first-passage flow
    mode2.S = sum( diag(diag(Q) - 1) %*% (ginv(diag(diag(Q))) %*% Q) %*%  diag(z) * dt) # storage from compartment-wise dissipative flow
    mode3.S = mode1.S  # storage from first-passage flow (equal to mode 1 at Steady state)
    mode4.S = sum(unpack(x)$y * dt)  # storage from boundary output flow

    #' packing up network statistics for output
    ns <- cbind('TSS'=TSS,
                'CIS'=CIS,
                'BSI'= BSI,'DSI'= DSI,'ISI'= ISI,
                'ID.S'=ID.S, 'ID.S.I'=ID.S.I,'ID.S.O'=ID.S.O,
                'HMG.S.O'=HMG.S.O,'HMG.S.I'=HMG.S.I,
                'NAS'=NAS,'NASP'=NASP,
                mode0.S,mode1.S,mode2.S,mode3.S,mode4.S)

                                        #'lam1P'=abs(lam1P),'rhoP'=abs(rhoP),
                                        #'lam1PP'=abs(lam1PP),'rhoPP'=abs(rhoPP),'AGG.S'=AGG.S)
                                        #re-orientation
    orient <- get.orient()
    if (orient == 'rc'){
      C <- t(C)
      P <- t(P)
      S <- t(S)
      VS <- t(VS)
      Q <- t(Q)
      CP <- t(CP)
      PP <- t(PP)
      SP <- t(SP)
      VSP <- t(VSP)
      QP <- t(QP)
    }else{}

    out <- list('X'=stor,'C'=C,'P'=P,'S'=S, 'VS'=VS,
                'Q'=Q,'CP'=CP,'PP'=PP,'SP'=SP, 'VSP'=VSP,
                'QP'=QP,'dt'=dt,'ns'=ns)

    return(out)
  }
}
#' enaStructure --- performes strucutral analysis of the
#' network graph (see Borrett et al. 2007)
#' INPUT = network object
#' OUTPUT = list of structure statistics
#'
#' S. Borrett and M. Lau | March 2011
#' ---------------------------------------------------







#' enaStructure --- performes strucutral analysis of the network graph (see
#' Borrett et al. 2007) INPUT = network object OUTPUT = list of structure
#' statistics
#' 
#' S. Borrett and M. Lau | March 2011
#' --------------------------------------------------- enaStructure ---
#' performes strucutral analysis of the network graph (see Borrett et al. 2007)
#' INPUT = network object OUTPUT = list of structure statistics
#' 
#' S. Borrett and M. Lau | March 2011
#' --------------------------------------------------- Structure Analyses of
#' Ecological Network
#' 
#' Analysis of the structure of an ecological flow network.
#' 
#' @param x A network object.
#' @return \item{A}{
#' 
#' } \item{ns}{A vector of structure based network statistics. These include n
#' = number of nodes, L = number of edges, C = connectivity, LD = link density,
#' ppr = pathway proliferation rate, lam1A = dominant eigenvalue, mlam1A =
#' multiplicity of dominant eigenvalue, rho = damping ratio, R = distance of
#' the dominant eigen value from the eigen spectra, d = difference between
#' dominant eigen value and link density, no.scc = number of strongly connected
#' components, no.scc.big = number of strongly connected components with more
#' than one node, pscc = percent of nodes in strongly connected components.  }
#' @author Matthew K. Lau Stuart R. Borrett
#' @seealso \code{\link{structure.statistics}}
#' @references Fath, B. D., Borrett, S. R. 2006. A Matlab function for Network
#' Environ Analysis.  Environ. Model. Softw. 21, 375-405.
#' @examples
#' 
#' 
#' 
#' data(troModels)
#' enaStructure(troModels[[6]])
#' 
#' 
#' 
#' @export enaStructure
enaStructure <- function(x = 'network object'){
                                        #Check for network class
  if (class(x) != 'network'){warning('x is not a network class object')}
  Flow <- t(as.matrix(x,attrname = 'flow')) #get flows
  A <- sign(Flow)   # get adjacency matrix
  sp <- structure.statistics(A)    # calls structure.statistics helper function
                                          #Output orientation
  orient <- get.orient()
  if (orient=='rc'){A <- t(A)}else{}
  return(list('A'=A,'ns'=sp))  # "A" is the adjacency matrix oriented
                                        # column to row and "sp" is a list of
                                        # structural network staistics
}
#' Trophic Aggregations (TroAgg) Analysis
#' 
#' It returns the data quantifying the underlying trophic structure of a given
#' model based on the interaction of the living and non-living nodes. It is
#' based on the Trophic Aggregations suggested by Lindeman (1942) and follows
#' the algorithm by Ulanowicz and Kemp (1979) implemented in NETWRK 4.2b. It
#' removes the Feeding cycles in the network beforehand to provide accurate
#' results.
#' 
#' 
#' @param x a network object.  This includes all weighted flows into and out of
#' each node. It should include separate respiration and export values for the
#' Canonical Exports and Canonical Respirations results respectively. It must
#' also include the "Living" vector that identifies the living (TRUE/FALSE)
#' status of each node. It must contain the non-living nodes at the end of the
#' node vector, the function \code{\link{netOrder}} can be used for the same.
#' @return \item{Feeding_Cycles}{List that gives the details of the Feeding
#' Cycles in the network. The output being according to the enaCycle function
#' applied to the Living components in the network} \item{A}{matrix that
#' distributes the species in integer Trophic Levels (Lindeman Transformation
#' Matrix). The dimension of A is (NL X NL) where NL is the number of Living
#' nodes.} \item{ETL}{vector of the Effective Trophic Level of each species.}
#' \item{M.flow}{vector of the Migratory flows, if present, in the network.}
#' \item{CI}{vector of Canonical Inputs to the integer trophic levels.
#' Displayed if the Migratory flows are present.} \item{CE}{vector of Canonical
#' exports or the exports from the integer trophic levels} \item{CR}{vector of
#' the Canonical Respirations or the respiration values for integer trophic
#' levels. } \item{GC}{vector of the input flow to a trophic level from the
#' preceeding trophic level. It represents the Grazing Chain for the network.}
#' \item{RDP}{vector of the Returns to Detrital Pool from each trophic level. }
#' \item{LS}{vector of the Lindeman trophic spine. It combines the Detrital
#' pool with the autotrophs and forms a monotonically decreasing sequence of
#' flows from one trophic level to the next, starting with the said
#' combination.} \item{TE}{vector of the trophic efficiencies i.e. the ratio of
#' input to a trophic level to the amount of flow that is passed on the next
#' level from it. } \item{ns}{vector of trophic aggregations based network
#' statistics. These include "Detritivory" the flow from the detrital pool to
#' the second trophic level, "DetritalInput" the exogenous inputs to the
#' detrital pool, "DetritalCirc" the circulation within the detrital pool,
#' "NCYCS" the number of feeding cycles removed, "NNEX" the number of feeding
#' cycle Nexuses removed and "CI" the Cycling Index for the Feeding Cycles.  }
#' @note This and other Ulanowicz school functions require that export and
#' respiration components of output be separately quantified.
#' 
#' This analysis involves the ENA Cycle analysis for removal of the Feeding
#' Cycles in the network. These are cycles amongst only the living nodes and
#' cause error in the trophic aggregations.
#' 
#' The analysis requires all the non-living nodes to be placed at the end in
#' the network object.
#' @author Pawandeep Singh
#' @seealso \code{\link{enaCycle}, \link{netOrder}}
#' @references %% ~put references to the literature/web site here ~ Lindeman,
#' R.L. 1942. The trophic-dynamic aspect of ecology. Ecology 23:399--418.
#' 
#' Ulanowicz, R.E. and Kemp, W.M.  1979. Towards canonical trophic
#' aggregations. The American Naturalist. 114:871--883.
#' 
#' Ulanowicz, R.E. 1995. Ecosystem trophic foundations: Lindeman exonerata. pp.
#' 549--560. B.C. Patten and S.E. Jorgensen (eds.) Complex Ecology: The
#' part-whole relation in ecosystems. Prentice Hall, New Jersey.
#' 
#' Ulanowicz, R.E. and Kay, J.J. 1991. A package for the analysis of ecosystem
#' flow networks. Environmental Software 6:131 -- 142.
#' @examples
#' 
#' 
#' 
#' data(troModels)
#' tro6 <- enaTroAgg(troModels[[6]])
#' attributes(tro6)
#' 
#' 
#' 
#' @export enaTroAgg
enaTroAgg <- function (x){
  if (class(x) != "network") {
    stop("x is not a network class object")
  }

                                        # Initials

  liv <- x %v% "living"     ##Living vector
  nl = sum(liv)             ##No. of living nodes
  N <- length(liv)
                                        #Living vector check
  liv2<-rep(FALSE,N)
  liv2[1:nl]<-rep(TRUE,nl)
  if(identical(liv,liv2) == FALSE) {
      stop('Non-living nodes must be at the end of the list.')
  }

  flow <- as.matrix(x, attrname = 'flow')
  XCHNGE<-flow
  Feeding_Cycles   <- cycliv(x)
  XCHNGE[1:nl,1:nl] <- Feeding_Cycles$ResidualFlows
  Ti <- x %v% "input"       ##AINPUT
  T <- Ti + apply(XCHNGE, 2, sum)

  exp<-x %v% "export"; exp[is.na(exp)] <- 0
  res<-x %v% "respiration"; res[is.na(res)] <- 0
  # ---------------------------------------------------------

                                        # Determining In-Migration and Obligate Producers
  BINPUT <- x %v% 'input'
  if(identical(XCHNGE,flow[1:nl,1:nl])) {print('Cycle free feeding transfers')}
  NMIG <- NPRM <- 0 ###NMIG - In-migration of Heterotrophs; NPRM - no. of obligate primary producers
  CANON <- rep(0,N)
  ###A compartment is an obligate producer iff it receives sustenance from no other compartment
  for(NP in 1:nl) {
  	if(Ti[NP]<=0){next}
  	for(i in 1:N) {
            ##Otherwise the input represents in-migration
  		if(XCHNGE[i,NP]<=0){next}
  		NMIG=NMIG+1
                ##Record the Migratory Input Temporary in CANON for use below
  		CANON[NP]=Ti[NP]
  		break
  	}
  	NPRM = NPRM+1
  }
  if(NPRM<=0){
  	warning("No Unambiguous primary producers found! All inputs assumed to be primary production!")
  	NMIG <- 0
  	CANON<- rep(0,N)
  	}
  mig.input<-rep(0,nl)
  if(NMIG>0) {
  	###Migratory Inputs. To be treated as Non-Primary Inflows
  	mig.input <- CANON[1:nl]
  }
  #---------------------------------------------------------

                                        #Recreate the matrix of feeding coefficients without the migratory inputs
  TL <- rep(0,N)
  FEED <- flow*0
  for(i in 1:N) {
  	##Store Throughputs in Vector TL
  	TL[i]<-0
  	for(j in 1:N) {TL[i]<-TL[i]+XCHNGE[j,i]}
  	TL[i]<-TL[i]+Ti[i]-CANON[i]
  	if(TL[i]<=0){TL[i]=1}
  	for(j in 1:N) {FEED[j,i]<-XCHNGE[j,i]/TL[i]}
  	BINPUT[i]<-(Ti[i]-CANON[i])/TL[i]
  }
  #---------------------------------------------------------

                                        #Create Lindeman Transformation Matrix
  CANON<- rep(0,N)
  A<-flow*0
  A[1,1:nl]<-BINPUT[1:nl]
  for(k in 2:nl) {
  	KM1=k-1
  	for(l in 1:nl) {
  		if((k<=2)&&(nl<N)) {
  			for(K2 in (nl+1):N) {A[k,l]<-A[k,l]+FEED[K2,l]}
  		}
  		for(j in 1:nl) {
  			A[k,l]<-A[k,l]+(A[KM1,j]*FEED[j,l])
  		}
  	}
  }
  if(nl<N) {for (i in (nl+1):N) {A[N,i]<-1}}
  ## A is the required Lindeman transformation matrix
  rownames(A) <- 1:N

                                        # 2. Effective Trophic Levels
  MF = matrix(1:N, nrow=N, ncol=N, byrow = 'FALSE')
  etl = rep(1,N)
  etl[1:nl] = apply((MF*A)[1:nl,1:nl],2,sum)


  ci <- Ti
  ci <- A %*% Ti
  ci <- as.vector(ci)


                                        # 3. Canonical Exports
  cel = exp
  ce = A %*% exp
  ce1 = as.vector(ce)

                                        # 4. Canonical Respirations
  crl = res
  cr = A%*%res
  cr1=as.vector(cr)

                                        # 5. Grazing Chain
  gc <- rep(0,nl)
  gc[1] <- sum(A[1,]*T)
  AT = A %*% flow %*% t(A)
  gc[2:nl]=apply(AT[1:(nl-1),1:nl,drop=FALSE],1,sum)

                                        # 6. Returns to Detrital Pool
  rtd <- AT[1:nl,N]
  rtd <- as.vector(rtd)

                                        # 7. Detrivory
  dtry <- sum(AT[N,1:nl])

                                        # 8. Input to Detrital Pool
  U <- A %*% Ti
  dinp<-0
  if(nl<N) {  dinp <- sum(U[(nl+1):N]) }

                                        # 9. Circulation within Detrital Pool
  dcir <- AT[N,N]

                                        # 10. Lindeman Spine
  ls = gc

  ls[1] = sum(rtd[2:nl]) + gc[1] + dinp
  ls[2] = gc[2]+dtry
                                        # 11. Trophic Efficiencies
  te=ls
  for(i in 1:nl){
    if(te[i]<=0){break}
    te[i]=te[i+1]/te[i]
  }
  te[is.na(te)] <- 0

                                        # Output Listing
  Detrivory<-dtry; DetritalInput<-dinp; DetritalCirc<-dcir
  ns <- cbind(Detrivory, DetritalInput, DetritalCirc, Feeding_Cycles$ns)
  if(NMIG>0) {
  	out <- list(Feeding_Cycles=Feeding_Cycles[1:(length(Feeding_Cycles)-1)], A = A[1:nl,1:nl], ETL = etl, M.Flow = mig.input, CI = ci, CE = ce1, CR = cr1, GC = gc, RDP = rtd, LS = ls,TE = te, ns=ns)
  }
  else{
  	out <- list(Feeding_Cycles=Feeding_Cycles[1:(length(Feeding_Cycles)-1)], A = A[1:nl,1:nl], ETL = etl, CE = ce1, CR = cr1, GC = gc, RDP = rtd, LS = ls,TE = te, ns=ns)

  }

  return(out)

  }#End of Function troAgg




#' enautility --- utility analysis of a flow network
#' INPUT = network object
#' OUTPUT = list of utility statistics
#'
#' M. Lau | July 2011
#' ---------------------------------------------------







#' enautility --- utility analysis of a flow network INPUT = network object
#' OUTPUT = list of utility statistics
#' 
#' M. Lau | July 2011 ---------------------------------------------------
#' enautility --- utility analysis of a flow network INPUT = network object
#' OUTPUT = list of utility statistics
#' 
#' M. Lau | July 2011 ---------------------------------------------------
#' Utility Analysis of Ecological Networks
#' 
#' Performs the flow and storage based utility analysis developed for
#' input-output network models of ecosystems.  It returns a set of matrices for
#' the direct and integral utilities as well as a set of utility based network
#' statistics.
#' 
#' @param x a network object.  This includes all weighted flows into and out of
#' each node.  For the storage utility analysis this must also include the
#' amount of energy--matter stored at each node (biomass).
#' @param type Determines whether the flow or storage utility analysis is
#' returned.
#' @param eigen.check LOGICAL: should the dominant eigenvalue be checked.  Like
#' enaFlow and enaStorage analyses, enaUtility analysis considers the utility
#' propigated over path lengths ranging for zero to infinity.  For utility
#' analysis to work properly, the path sequence must converge.  enaUtility
#' checks to see if the utility path sequence is convergent by finding the
#' dominant eigenvalue of the direct utility matrix.  If this eigenvalue is
#' less than 1, the sequence is convergent and the analysis can be applied; if
#' the dominant eigenvalue is greater than one, then the anlysis cannot be
#' applied.  By default, the function will not return utility values if the
#' eigenvalue is larger than one; however, if eigen.check is set to FALSE, then
#' the function will be applied regardless of the mathematic validity.
#' @param balance.override LOGICAL: should model balancing be ignored.
#' enaUtility assumes that the network model is at steady-state.  The default
#' setting will not allow the function to be applied to models not at
#' steady-state.  However, when balance.override is set to TRUE, then the
#' function will work regardless.
#' @param tol The integral utility matrix is rounded to the number of digits
#' specified in tol.  This approximation eleminates very small numbers
#' introduced due to numerical error in the ginv function.  It does not
#' eliminate the small numerical error introduced in larger values, but does
#' truncate the numbers.
#' @return \item{D}{Direct flow utility intensity matrix.  (fij-fji)/Ti for
#' i,j=1:n} \item{U}{Nondimensional integral flow utility} \item{Y}{Dimensional
#' integral flow utility} \item{ns}{If type is set to 'flow', this is a list of
#' flow utility network statistics including: the dominant eigenvalue of D
#' (lambda\_1D), flow based network synergism (synergism.F), and flow based
#' network mutualism (mutualism.F).} \item{DS}{Direct storage utility intensity
#' matrix.  (fij-fji)/xi for i,j=1:n} \item{US}{Nondimensional integral storage
#' utility} \item{YS}{Dimensional integral storage utility} \item{ns}{If type
#' is set to 'storage', this is a list of storage utility network statistics
#' including: the dominant eigenvalue of DS (lambda_1DS), storage based network
#' synergism (synergism.S), and storage based network mutualism (mutualism.S).}
#' @author Matthew K. Lau Stuart R. Borrett
#' @seealso \code{\link{enaFlow},\link{enaStorage},\link{enaMTI}}
#' @references Fath, B.D. and Patten, B.C. 1998. Network synergism: emergence
#' of positive relations in ecological systems.  Ecol. Model. 107:127--143.
#' 
#' Fath, B.D. and Borrett, S.R. 2006. A Matlab function for Network Environ
#' Analysis. Environ. Model. Soft. 21: 375--405.
#' 
#' Patten, B.C. 1991.  Network ecology: Indirect determination of the
#' life-environment relationship in ecosystems.  In: Higashi, M. and Burns, T.
#' (eds). Theoretical Studies of Ecosystems: The Network Perspective. Cambridge
#' University Press.  New York.
#' @examples
#' 
#' 
#' 
#' data(troModels)
#' U <- enaUtility(troModels[[6]], type = "flow", eigen.check = FALSE)
#' attributes(U)
#' US <- enaUtility(troModels[[6]], type = "storage", eigen.check = FALSE)
#' 
#' 
#' 
#' @export enaUtility
enaUtility <- function(x, type=c('flow','storage'),
                       eigen.check=TRUE,
                       balance.override=FALSE,tol=10){
                                        #Missing Data Check
    if (type == 'storage' && any(is.na(x%v%'storage'))){
        warning('This function requires quantified storage values.')
    }else{
        orient <- get.orient()
                                        #Check for network class
        if (class(x) != 'network'){warning('x is not a network class object')}

                                        #set default for type == 'flow'
        if (length(type) > 1){type <- 'flow'}

                                        #Check for balancing
        if (balance.override){}else{
            if (any(list.network.attributes(x) == 'balanced') == FALSE){x%n%'balanced' <- ssCheck(x)}
            if (x%n%'balanced' == FALSE){warning('Model is not balanced'); stop}
        }
                                        #unpack data from x
        Flow <- t(as.matrix(x, attrname = 'flow')) #flows
        input <- x%v%'input' #inputs
        stor <- x%v%'storage' #storage values
        I <- Flow*0; diag(I) <- 1 #create the identity matrix
        T. <- input + apply(Flow,1,sum)
        FD <- Flow;
        diag(FD) <- -T.
                                        #
        if (type == 'flow'){
                                        #flow utilities
            D <- ginv(diag(T.)) %*% (FD - t(FD))
            rownames(D) <- colnames(D) <- colnames(Flow)

            if (eigen.check & abs(eigen(D)$values[1]) > 1){
                print(paste('Largest eigen value of D > 1:',eigen(D)$values[1]),quote=FALSE)
                out <- NA
            }else{
                U <- ginv(I-D) #non-dimensional integral flow utility
                U <- round(U,tol)
                Y <- diag(T.) %*% U #dimensional integral flow utility
                rownames(U) <- colnames(U) <- colnames(Flow)
                rownames(Y) <- colnames(Y) <- colnames(Flow)

                                        #indices
                synergism.F <- bcratio(Y) #flow benefit cost ratio (calls other function) (Synergism)
                mutualism.F <- bcratio(sign(Y)) # flow ratio of positive to negative signs )

                # get relational data
                R <- relationalChange(D,Y)
                SD <- R$Direct.Signs
                SY <- R$Integral.Signs
                R.table <- R$Relations.Table
                names(R.table) <- c("From","To","Direct","Integral","changed")
                rownames(R.table) <- c(1:dim(R.table)[1])
                                        #re-orient
                if (orient == 'rc'){
                    D <- t(D)
                    U <- t(U)
                    Y <- t(Y)
                    SD <- t(SD)
                    SY <- t(SY)
                }else{}

                ns <- cbind('lam1D' = abs(eigen(D)$values[1]),
                            'relation.change.F' = R$ns[1],
                            'synergism.F' = synergism.F,
                            'mutualism.F' = mutualism.F)

                out <- list('D'=D, 'SD' = SD,
                            'U'=U,
                            'Y'=Y, 'SY' = SY,
                            'Relations.Table' = R.table,
                            'ns'=ns) #pack output

            }

            # ------------------------------------------------------------------------
        }else if (type == 'storage'){
                                        #storage utilities
            x <- stor
            DS <- ginv(diag(x)) %*% (FD - t(FD))
            rownames(DS) <- colnames(DS) <- colnames(Flow)

            if (eigen.check & abs(eigen(DS)$values[1]) > 1){
                print(paste('Largest eigen value of DS > 1:',eigen(DS)$values[1]),quote=FALSE)
                out <- NA
      }else{

          US <- ginv(I - DS)
          US <- round(US,tol)
          YS <- diag(T.) %*% US
          rownames(US) <- colnames(US) <- colnames(Flow)
          rownames(YS) <- colnames(YS) <- colnames(Flow)


                                        #indices
          synergism.S <- bcratio(YS) #storage benefit cost ratio (calls other function) (Synergism)
          mutualism.S <- bcratio(sign(YS)) #storage ratio of positive to negative signs (Y/abs(Y) == sign of Y)

                          # get relational data
          R <- relationalChange(DS,YS)
          SD <- R$Direct.Signs
          SY <- R$Integral.Signs
          R.table <- R$Relations.Table
          names(R.table) <- c("From","To","Direct","Integral","changed")
          rownames(R.table) <- c(1:dim(R.table)[1])

                                        #re-orient
          if (orient == 'rc'){
              DS <- t(DS)
              US <- t(US)
              YS <- t(YS)
              SD <- t(DS)
              SY <- t(SY)

          }else{}

          ns <- cbind('lam1DS'=abs(eigen(DS)$values[1]),
                      'relation.change.S' = R$ns[1],
                      'synergism.S' = synergism.S,
                      'mutualism.S'=mutualism.S)

          out <- list('DS'=DS,'SD'=SD,
                      'US'=US,
                      'YS'=YS,'SY'=SY,
                      'Relations.Table' = R.table,
                      'ns'=ns) #package output
      }
        }
                                        #labeling
        if (length(out)>1){
            for (i in 1:(length(out)-1)){
                if (class(out[[i]])=='matrix'){
                    rownames(out[[i]]) <- colnames(out[[i]]) <- colnames(Flow)
                }
            }
        }
                                        #output
        return(out)
    }
}
#' environCentrality --- calculates the centrality of 
#' flow network environs
#' INPUT = environ matrix
#' OUTPUT = in-going, out-going and average centralities
#' 
#' M. Lau | July 2011
#' ---------------------------------------------------







#' environCentrality --- calculates the centrality of flow network environs
#' INPUT = environ matrix OUTPUT = in-going, out-going and average centralities
#' 
#' M. Lau | July 2011 ---------------------------------------------------
#' environCentrality --- calculates the centrality of flow network environs
#' INPUT = environ matrix OUTPUT = in-going, out-going and average centralities
#' 
#' M. Lau | July 2011 ---------------------------------------------------
#' Environ Centrality an Ecological Network
#' 
#' This function calculates the input, output, and average environ centrality
#' of the nodes in the network (Fath and Borret, 2012).  This is a type of
#' weighted degree centrality that indicates the relative importance of the
#' nodes in the flow activity in the network.
#' 
#' @param x A square matrix.  Usually the integral flow marix from enaFlow. The
#' assumption is that the flows are oriented column to row.
#' @return \item{ECin}{input oriented environ centrality} \item{ECout}{output
#' oriented environ centraility} \item{AEC}{average environ centrality (average
#' of input and output)}
#' @author Matthew K. Lau Stuart R. Borrett
#' @seealso \code{\link{enaFlow}}
#' @references Fann, S.L. and Borrett, S.R. 2012. Environ centrality reveals
#' the tendency of indirect effects to homogenize the functional importance of
#' species in ecosystems.  Journal of Theoretical Biology 294: 74-86.
#' @examples
#' 
#' 
#' 
#' data(troModels)
#' F <- enaFlow(troModels[[6]])
#' ec <- environCentrality(F$N)
#' attributes(ec)
#' barplot(sort(ec$AEC, decreasing = TRUE), col = 4, ylab = "Average Environ Centrality", 
#'     ylim = c(0, 0.4))
#' 
#' 
#' 
#' @export environCentrality
environCentrality <- function(x='matrix'){
  if (class(x) != 'matrix'){warning('x is not a matrix class object')}
  ECin <- rowSums(x)/sum(rowSums(x))
  ECout <- colSums(x)/sum(rowSums(x))
  AEC <- (ECin + ECout)/2
  names(ECin) <- rownames(x)
  names(ECout) <- rownames(x)
  names(AEC) <- rownames(x)
  return(list('ECin'=ECin,'ECout'=ECout,'AEC'=AEC))
}
#' findPathLength --- calculates the flows over a 
#' sequence up to a maximum path length
#' INPUT = network object
#' OUTPUT = a list of flow statistics over paths
#' 
#' S. Borrett and M. Lau | July 2011
#' ---------------------------------------------------







#' findPathLength --- calculates the flows over a sequence up to a maximum path
#' length INPUT = network object OUTPUT = a list of flow statistics over paths
#' 
#' S. Borrett and M. Lau | July 2011
#' --------------------------------------------------- findPathLength ---
#' calculates the flows over a sequence up to a maximum path length INPUT =
#' network object OUTPUT = a list of flow statistics over paths
#' 
#' S. Borrett and M. Lau | July 2011
#' --------------------------------------------------- Cumulative Flow over a
#' Range of Path Lengths
#' 
#' Calculates the flow throughout the entire network over a given path length.
#' 
#' @param x Network model object.
#' @param maxPath The maximum path length to calculate total flow.
#' @param plot.sw LOGICAL: should a plot be generated showing flow
#' accumulation?
#' @return \item{thresholds}{thresholds indicating the development of
#' throughflow as path length increases: the path length at which indirect flow
#' exceeds direct flow (mID), path length at which 50\%, 90\%, and 95\% of
#' total system throughflow is achieved (m50, m90, and m95, respectively)}
#' \item{tf}{total flow across paths from length 0 (Boundary inputs) to
#' maxPath} \item{ctf}{cumulative total flow from path length 0 to maxPath}
#' @author Matthew K. Lau Stuart R. Borrett
#' @seealso \code{\link{enaFlow}}
#' @references Borrett, S.R, Patten, B.C., Whipple, S.J. 2010.  Rapid
#' development of indirect effects in ecological networks.  Oikos
#' 119:1136--1148.
#' @examples
#' 
#' 
#' 
#' data(troModels)
#' pl10 <- findPathLength(troModels[[6]], plot.sw = TRUE, maxPath = 10)
#' names(pl10)
#' pl10$thresholds
#' 
#' 
#' 
#' @export findPathLength
findPathLength <- function(x,maxPath=100,plot.sw=FALSE){
  ##
  if(ssCheck(x)=="FALSE"){x = balance(x)}  # ensure the models is balanced
  oo <- get.orient() #original orientation
  if (oo == 'school'){oo <- 'internal'}
  set.orient('internal')
  Flow <- enaFlow(x)   # perform flow analysis
  set.orient(oo)
                                        #
  TST <- Flow$ns[2]
                                        # find Total Flow over each path length
  k <- 0:maxPath
  tf <- unlist(lapply(k,function(k) sum( mExp(Flow$G,k) %*% as.matrix(x%v%'input'))))
  tfi <- tf/TST # total flow intensity flow/TST
                                        # find cumulative flow percentage
  k <- 1:(maxPath+1)
  ctf <- unlist(lapply(k, function(k) sum(tfi[1:k])))
                                        # find thresholds
  m50 <- (min(which(ctf>=0.5))-1) # need to subtract 1 becuase index 1 is path length 0.
  m90 <- (min(which(ctf>=0.9))-1)
  m95 <- min(which(ctf>=0.95))-1

  if(Flow$ns[8]>1){
                                        # find cumulative indirect flow
    direct <- tf[2]   # k =1 is boundary, k = 2 is direct
    k <- 3:(maxPath+1)
    cindirect <- unlist(lapply(k, function(k) sum(tf[3:k]))) 
    mID <- min(which(cindirect>direct))+1
  } else {mID <- NA}
    
  if(plot.sw){
    opar <- par(las=1)
    plot(0:(length(ctf)-1),ctf,type="b",pch=20,col="blue",ylim=c(0,1),
         xlab="Path Length",ylab="Cumulative Flow Intensity",axes=FALSE)
    axis(2,at=c(0,0.25,0.5,0.75,0.95,1))
    axis(1,at=c(seq(0,maxPath,by=maxPath/5),m50,m95))
    box()
    points(c(-m50,m50),c(0.5,0.5),type="l",lty=2)
    points(c(m50,m50),c(-m50,0.5),type="l",lty=2)
    points(c(-m95,m95),c(0.95,0.95),type="l",lty=2)
    points(c(m95,m95),c(-m95,0.95),type="l",lty=2)
    par(opar)
    rm(opar)
  }
  thresholds <- c("mID"=mID,"m50"=m50,"m90"=m90,"m95"=m95)
 return(list("thresholds"=thresholds,"tf"=tf,"ctf"=ctf))
  
}
#' force.balance --- repeatedly applies balance until 
#' sub-tolerance is reached
#' INPUT = network model
#' OUTPUT = balanced model
#' M. Lau 1 Oct 2012
#' ---------------------------------------------------







#' force.balance --- repeatedly applies balance until sub-tolerance is reached
#' INPUT = network model OUTPUT = balanced model M. Lau 1 Oct 2012
#' --------------------------------------------------- force.balance ---
#' repeatedly applies balance until sub-tolerance is reached INPUT = network
#' model OUTPUT = balanced model M. Lau 1 Oct 2012
#' --------------------------------------------------- Repeated Application the
#' Balance Function
#' 
#' This function repeatedly balances a model, sequentially with the output
#' being passed back to the balance function, until it is within tolerance or
#' the maximum number of iterations is reached.
#' 
#' 
#' @param x A network object.
#' @param tol Percent error tolerance for difference between inputs and
#' outputs.
#' @param max.itr Maximum number iterations.
#' @param method The balancing method to use, see balance. DEFAULT = AVG2.
#' @return Returns a balanced network model.
#' @author Matthew K. Lau Stuart R. Borrett
#' @seealso \code{\link{balance}}
#' @references Allesina, S., Bondavalli, C., 2003.Steady state of ecosystem
#' flow networks: a comparison between balancing procedures.Ecological
#' Modelling 165(2-3):231-239.
#' @examples
#' 
#' 
#' 
#' data(troModels)
#' ssCheck(troModels[[1]])
#' fb.model = force.balance(troModels[[2]])  #produces a balanced model
#' 
#' 
#' 
#' @export force.balance
force.balance <- function(x,tol=5,max.itr=10,method='AVG2'){
  n.itr <- 1 # initiate counter
  while(ssCheck(x)==FALSE & n.itr<max.itr){
    x <- balance(x,method=method,tol=tol)
    n.itr <- n.itr + 1
  }
  if (n.itr>=max.itr){
    warning('Maximum iterations reached.')
  }else{
    return(x)
  }
}
#' get.ns.R
#' Input = network model
#' Output = a vector of global network statistics from ena
#'
#' Borrett | July 4, 2012
#' -----------------------------------







#' get.ns.R Input = network model Output = a vector of global network
#' statistics from ena
#' 
#' Borrett | July 4, 2012 ----------------------------------- get.ns.R Input =
#' network model Output = a vector of global network statistics from ena
#' 
#' Borrett | July 4, 2012 ----------------------------------- Quick Calculation
#' of a Range of Network Statistics.
#' 
#' This is a high level function for calculated the main network analyses
#' (Ascendancy, Flow, Structure, Storage and Utility) on an ecological network.
#' 
#' @param x A network object.
#' @param balance.override Turns off balancing and balance checking.
#' @return Returns the network statistics (ns) of all of the major ENA
#' functions: enaStructure, enaFlow, enaAscendency, enaStorage and enaUtility
#' (both flow and storage).
#' @author Matthew K. Lau Stuart R. Borrett David E. Hines
#' @seealso
#' \code{\link{enaStructure}},\code{\link{enaFlow}},\code{\link{enaAscendency}},\code{\link{enaUtility}}
#' @references Fath, B. D., Borrett, S. R. 2006. A Matlab function for Network
#' Environ Analysis.  Environ. Model. Softw. 21, 375-405.
#' @examples
#' 
#' 
#' 
#' data(troModels)
#' get.ns(troModels[[6]])
#' 
#' 
#' 
#' @export get.ns
get.ns <- function(x,balance.override=FALSE){
                                        #Check for network class
  if (class(x) != 'network'){warning('x is not a network class object')}

                                        #Check for balancing
  if (balance.override){}else{
    if (any(list.network.attributes(x) == 'balanced') == FALSE){x%n%'balanced' <- ssCheck(x)}
    if (x%n%'balanced' == FALSE){warning('Model is not balanced'); stop}
  }
  # runs selected ena analyses that return global network statistics
  st <- enaStructure(x)$ns
  Flow <- enaFlow(x)$ns
#  asc <- enaAscendency(x)  # enaFlow now includes the ascendency measures
  s <- enaStorage(x)$ns
  u.f <- enaUtility(x,type='flow',eigen.check=FALSE)$ns
  u.s <- enaUtility(x,type='storage',eigen.check=FALSE)$ns
  ns <- data.frame(st,Flow,s,u.f,u.s)
  rownames(ns) <- ""
  return(ns)
}

#' get.orient --- returns the global orientation
#' INPUT = none
#' OUTPUT = returns the current orientation of matrices
#' 
#' M. Lau | 16 Jun 2013
#' ---------------------------------------------------







#' get.orient --- returns the global orientation INPUT = none OUTPUT = returns
#' the current orientation of matrices
#' 
#' M. Lau | 16 Jun 2013 ---------------------------------------------------
#' get.orient --- returns the global orientation INPUT = none OUTPUT = returns
#' the current orientation of matrices
#' 
#' M. Lau | 16 Jun 2013 --------------------------------------------------- Get
#' the Current Global Matrix Orientation Setting
#' 
#' Returns the current setting for the expected orientation of all matrices,
#' which is either 'rc' (DEFAULT) or 'school' (output orientation as expected
#' for the school of analysis for a given function).
#' 
#' This function is intended to provide increase flexibility for users of both
#' the Patten and Ulanowicz schools of ENA.
#' 
#' @author M.K. Lau and S.R. Borrett
#' @export get.orient
get.orient <- function(){
  current.orientation <- get('orientation',envir=environment(set.orient))
  return(current.orientation)
}
#'# mExp --- calculate the exponent of a given matrix
#'# INPUT = a matrix (x) and the exponent (n)
#'# OUTPUT = the resulting exponentiated matrix
#'# 
#'# Alberto Monteiro (https://stat.ethz.ch/pipermail/
#'# r-help/2007-May/131330.html)
#'# ___________________________________________________







#' # mExp --- calculate the exponent of a given matrix # INPUT = a matrix (x)
#' and the exponent (n) # OUTPUT = the resulting exponentiated matrix # #
#' Alberto Monteiro (https://stat.ethz.ch/pipermail/ #
#' r-help/2007-May/131330.html) #
#' ___________________________________________________ # mExp --- calculate the
#' exponent of a given matrix # INPUT = a matrix (x) and the exponent (n) #
#' OUTPUT = the resulting exponentiated matrix # # Alberto Monteiro
#' (https://stat.ethz.ch/pipermail/ # r-help/2007-May/131330.html) #
#' ___________________________________________________ Calculates the Exponent
#' of a Matrix
#' 
#' Function for calculating the pathway proliferation of flows in a network
#' model through matrix exponentiation.
#' 
#' 
#' @param x A matrix.
#' @param n Desired exponent (i.e. the path length).
#' @return Returns an exponentiated flow matrix.
#' @author Alberto Monteiro
#' (https://stat.ethz.ch/pipermail/r-help/2007-May/131330.html) Matthew K. Lau
#' @seealso \code{\link{findPathLength}}
#' @references This function was originally designed by Alberto Monteiro in the
#' following R help thread:
#' https://stat.ethz.ch/pipermail/r-help/2007-May/131330.html.
#' @export mExp
mExp <- function(x='matrix', n=2){
  if (n == 1) return(x)
  result <- diag(1, ncol(x))
  while (n > 0) {
    if (n %% 2 != 0) {
      result <- result %*% x
      n <- n - 1
    }
    x <- x %*% x
    n <- n / 2
  }

  rownames(result) <- colnames(result)

  return(result)

}
#'## Function to order the nodes in a Network in enaR
#'## Singh P. | July 2014
#'## -----------------------------------------







#' ## Function to order the nodes in a Network in enaR ## Singh P. | July 2014
#' ## ----------------------------------------- ## Function to order the nodes
#' in a Network in enaR ## Singh P. | July 2014 ##
#' ----------------------------------------- Reorder Nodes in a Network in enaR
#' 
#' Reorders nodes in a network either through a user defined node order vector
#' or by default places the non-living nodes to the end of the node vector,
#' minimizing the order change for other nodes.
#' 
#' 
#' @param x A network object. This includes all weighted flows into and out of
#' each node.
#' @param order An integer vector of length N, where N is number of nodes in x,
#' specifying the new order of the nodes (by default order = 0, which indicates
#' moving non-living nodes to the end)
#' @return Returns a network object with nodes ordered as per the node order
#' vector or without the node order vector, by default moves the non-living
#' nodes to the end of the node vector, minimizing the order change for other
#' nodes.
#' @note The node order vector "order" must be of length equal to the number of
#' nodes in x (i.e. N) and must contain all integers from 1 to N.
#' 
#' This function can be used with default conditions (i.e. without "order"
#' vector) to reorder the nodes of a network which does not have non-living
#' nodes placed at the end so that the Trophic Aggregations analysis
#' (enaTroAgg) can be run on the reordered model.
#' @author Pawandeep Singh
#' @seealso \code{\link{enaTroAgg}}
#' @examples
#' 
#' 
#' 
#' data(troModels)
#' new.network <- netOrder(troModels[[6]], c(1, 3, 2, 5, 4))
#' # new.network is the required rearranged network with nodes in the desired order.
#' 
#' 
#' 
#' @export netOrder
netOrder <- function(x,order=0) {
    if (class(x) != "network") {
        stop("x is not a network class object")
    }
                                        # Load Initials
    flow <- as.matrix(x, attrname = "flow")
   # input <- x %v% "input"
   # resp <- x %v% "respiration"
   # export <- x %v% "export"
   # output <- x %v% "output"
   # storage <- x %v% "storage"
    living <- x %v% "living"
   # names <- x %v% "vertex.names"
    N <- length(living)
                                        # Determine Order (ordr)

    if(identical(order,0)==TRUE) {
        order<-rep(0,N)
        liv1<-which(living==TRUE)
        liv2<-which(living==FALSE)
        order<-c(liv1,liv2)
        if(identical(order,1:N)==TRUE) {warning('Network meets default conditions, no changes made')}
    }



                                        # Rearrange Network Characteristics
   # living <- living[ordr]
    flow  <- flow[order,order]
   # export <- export[ordr]
   # resp <- resp[ordr]
   # storage <- storage[ordr]
   # output <- output[ordr]
   # input <- input[ordr]
   # names <- names[ordr]


                                        # Modify Network
    x<-permute.vertexIDs(x,order)
    set.edge.attribute(x, 'flow', flow[flow>0])
    #x %v% "input" <- input
    #x %v% "respiration" <- resp
    #x %v% "export" <- export
    #x %v% "output" <- output
    #x %v% "storage" <- storage
    #x %v% "living" <- living
    #x %v% "vertex.names" <- names


                                        # Return the ordered network
    return(x)

}
# pack --- helper function for inputing flow
# network information into a network object
# INPUT = flow network model components
# OUTPUT = a network object
# M.Lau & S.R. Borrett | July 2014
# ------------------------------------







#' Compile Network Information into a Network Class
#' 
#' This function provides a flexible framework for importing flow network
#' information into a network class object for analyses.
#' 
#' 
#' @param flow The flow matrix.
#' @param input The inputs into the system.
#' @param respiration The quantities respired from the system.
#' @param export The exports from the system.
#' @param output The output (i.e. exports + respiration) from the system.
#' @param storage The quantities stored in compartments within the system.
#' @param living A logical vector indicating whether a node is either 'living'
#' (= TRUE) or 'dead' (=FALSE).
#' @return Returns a network object for the supplied model.
#' @author Matthew K. Lau Stuart R. Borrett
#' @seealso \code{\link{unpack}}
#' @export pack
pack <- function(flow,input=NA,respiration=NA,export=NA,output=NA,storage=NA,living=NA){
                                        #Warn if missing both
  if (all(is.na(respiration)) & all(is.na(export))){
    warning('Missing or NA resipiration/export values.')
  }else if (any(is.na(output) == FALSE) & 
            any(is.na(export) == FALSE)){
                respiration <- output - export
            }else if (any(is.na(output) == FALSE) & 
                      any(is.na(respiration) == FALSE)){
                export <- output - respiration
            }else if (all(is.na(output)) & any(c(is.na(respiration), is.na(export)) == FALSE)){
                export. <- export
                respiration. <- respiration
                export.[is.na(export)] <- 0
                respiration.[is.na(respiration)] <- 0
                output <- export. + respiration.
                output[(is.na(export) + is.na(respiration)) == 2] <- NA
            }
                                        #Add rownames
  if (length(rownames(flow))==0){rownames(flow) <- colnames(flow) <- as.character(1:nrow(flow))}
                                        #Compiling the objects into a list
  x <- list('flow' = as.matrix(flow),'input' = input,'export' = export,
            'respiration' = respiration, 'storage' = storage,'living'=living)

                                        #Warning for missing components
  if(any(is.na(unlist(x)))){
    missing <- print(names(unlist(x))[is.na(unlist(x))])
    if (length(missing)>1){
      for (i in 2:length(missing)){
        missing[1] <- paste(missing[1],missing[i],sep=', ')
      }
    }
    warning(paste('Missing model components:',missing[1],sep=' '))
  }else{}
                                        #initializing the network object using the flow matrix
  y <- network(x[[1]],directed=TRUE,loops=TRUE)
  # edge
  set.edge.attribute(y,names(x)[1],as.numeric(x[[1]]))
                                        #packing up the attributes into the network object (y)
  flow <- as.matrix(flow)
  rownames(flow) <- colnames(flow)
  set.edge.attribute(y,'flow',flow[flow>0])
  # vertex
  set.vertex.attribute(y,'input',input)
  set.vertex.attribute(y,'export',export)
  set.vertex.attribute(y,'respiration',respiration)
  set.vertex.attribute(y,'output',output)
  set.vertex.attribute(y,'storage',storage)
  set.vertex.attribute(y,'living',living)
  set.vertex.attribute(y,'vertex.names',rownames(flow))
                                        #naming the rows and columns of the flow matrix and storing
                                        #it in the network attributes



                                        #check if model is balanced
  y %n% 'balanced' <- ssCheck(y) #check if the model is balanced

  return(y)
}
### Reading econet models
### 1 Mar 2016
### mklau







#' Read an EcoNet model.
#' 
#' This function allows the user to access models that are formatted for
#' EcoNet, the web-based interface for conducting ENA
#' (http://eco.engr.uga.edu/), by Caner Kazanci at the University of Georgia.
#' 
#' 
#' @param x An object with the EcoNet formatted file, which can be read into R
#' using readLines.
#' @param verbose LOGICAL: should warnings be suppressed?
#' @return Returns the model formatted as a network object.
#' @author Matthew K. Lau
#' @seealso \code{\link{EcoNetWeb}}
#' @references Kazanci, C., 2007. EcoNet: A new software for ecological
#' modeling, simulation and network analysis, Ecol. Model., Vol 208/1 pp 3-8.
#' @export read.EcoNet
read.EcoNet <- function(x,verbose=FALSE){
    if (!(verbose)){options(warn=-1)}
    x <- x[!(grepl('<',x)) & grepl('=',x)]
    x <- x[!(grepl('\\#',x))]
    if (any(!(grepl('c=',x)) & grepl('=',x))){
        stor <- x[grepl('=',x) & !(grepl('c=',x))]
    }else{stor <- NA}
    x <- x[grepl('=',x) & grepl('c=',x)]
    x <- paste(x,collapse=';')
    x <- strsplit(x,' ')[[1]]
    x <- x[x != '']
    x <- paste(x,collapse='')
    x <- strsplit(x,';')[[1]]
    flo <- x[!(grepl('\\*',x))]
    inp <- x[(grepl('\\*->',x))]
    out <- x[(grepl('->\\*',x))]
    flo <- sub('->',';',flo)
    flo <- sub('c=',';',flo)
    flo <- do.call(rbind,strsplit(flo,split=';'))
    flow <- matrix(0,nrow=length(unique(c(flo[,1],flo[,2]))),
                   ncol=length(unique(c(flo[,1],flo[,2]))))
    rownames(flow) <- colnames(flow) <- unique(c(flo[,1],flo[,2]))
    for (i in 1:nrow(flo)){
        flow[rownames(flow) == flo[i,1],colnames(flow) == flo[i,2]] <- as.numeric(flo[i,3])
    }
    inp <- do.call(rbind,strsplit(sub('\\*->','',inp),'c='))
    input <- as.numeric(inp[,2])
    names(input) <- inp[,1]
    out <- do.call(rbind,strsplit(sub('->\\*','',out),'c='))
    output <- as.numeric(out[,2])
    names(output) <- out[,1]
    if (is.na(stor[[1]])){
        storage <- rep(0,nrow(flow))
    }else{
        stor <- paste(stor,collapse='')
        stor <- strsplit(stor,split='')[[1]]
        stor[!(stor %in% c(LETTERS,letters,0:9,'=','.','_',' '))] <- ','
        stor <- paste(stor,collapse='')
        stor <- paste(strsplit(stor,' ')[[1]],collapse='')
        stor <- strsplit(stor,',')[[1]]
        stor <- strsplit(stor,'=')
        stor <- do.call(rbind,stor)
        storage <- as.numeric(stor[,2])
        names(storage) <- stor[,1]
    }
    if (length(input) != nrow(flow)){
        inp <- rep(0,(nrow(flow) - length(input)))
        names(inp) <- rownames(flow)[!(rownames(flow) %in% names(input))]
        input <- c(input,inp)
    }
    if (length(output) != nrow(flow)){
        outp <- rep(0,(nrow(flow) - length(output)))
        names(outp) <- rownames(flow)[!(rownames(flow) %in% names(output))]
        output <- c(output,outp)
    }
    input <- input[match(names(input),rownames(flow))]
    output <- output[match(names(output),rownames(flow))]
    storage <- storage[match(names(storage),rownames(flow))]
    return(pack(flow=flow,input=input,output=output,storage=storage))
}
#' R function to read in a matrix formatted as Mdloti (Ursula Sharler)
#' Borrett | Sept. 12, 2012, MKL July 2013
#' Updated - Borrett, May 2016 - to use pack() to create the network data object.
#' ------------------------







#' R function to read in a matrix formatted as Mdloti (Ursula Sharler) Borrett
#' | Sept. 12, 2012, MKL July 2013 Updated - Borrett, May 2016 - to use pack()
#' to create the network data object. ------------------------ R function to
#' read in a matrix formatted as Mdloti (Ursula Sharler) Borrett | Sept. 12,
#' 2012, MKL July 2013 Updated - Borrett, May 2016 - to use pack() to create
#' the network data object. ------------------------ Read ENA Model from an
#' Mdloti Formatted Excel File
#' 
#' This function reads network data from an excel file commonly used by Ursula
#' Sharler.  The file has three header lines (name/source, number of
#' compartments, number of living nodes) and then a n+2 x n+2 matrix of flows.
#' This is the flow matrix with an additional row for imports and biomass each
#' and additional columns for exports and respirations.
#' 
#' 
#' @param file The name and path for the data file.  This function assumes the
#' data are stored on the first sheet of an Microsoft Excel formatted. NOTE:
#' this function depends on the read.xlsx function from the xlsx package, which
#' requires that the entire path be specified from the root directory (i.e. the
#' absolute path).
#' @return Returns the network object.
#' @author Stuart R. Borrett
#' @seealso \code{\link{read.scor}}
#' @references Fath, B. D., Borrett, S. R. 2006. A Matlab function for Network
#' Environ Analysis.  Environ. Model. Softw. 21, 375-405.
#' @export read.enam
read.enam<- function(file="file path and name"){
                                        #I have assumed the file is formatted as an excel speadsheet.
                                        #The data must be on the first sheet in the workbook.
  x <- as.matrix(read.xls(file,sheet=1,header=FALSE))
  mname <- as.character(x[1,1]); # Get Model ID
  n <- as.numeric(as.character(x[2,2])) # number of nodes
  liv <- as.numeric(as.character(x[3,2])) # number of nodes
  a <- n+6+1 # ending row of flows matrix -- assumes Flows start on row 6 and Imports and Biomasses are at the end
  b <- n+2+2 # ending column of flows matrix -- assumes exports and respirations are at the end
  m <- x[6:a,3:b] # Matrix of Flows
  m <- apply(m,2,as.numeric)
  rownames(m) <- colnames(m) <- as.character(x[6:a,2]) # node names
  m[is.na(m)] <- 0 # replace NAs with zeros
  Flow <- m[1:n,1:n] # flow matrix
  imports <- m[(n+1),1:n]
  biomass <- as.numeric(unlist(m[(n+2),1:n]))
  exports <- as.numeric(unlist(m[1:n,(n+1)]))
  respiration <- as.numeric(unlist(m[1:n,(n+2)]))
  LIV <-  c(rep(TRUE,liv),rep(FALSE,n-liv))
                                       # packing up the attributes into the network object (y)
  y <- pack(flow = Flow, input = imports, respiration = respiration, export = exports, storage = biomass, living = LIV)
  return(y)
}
#' read.nea.RData
#' INPUT = Model Data (flows, inputs, outputs, storage) formatted as for NEA.m, saved as CSV file
#'        S=  |[F][z][X]|
#'            |[y][0][0]|
#' OUPUT = R Network data object for use with enaR
#'
#' Borrett | July 15, 2013
#' --------------------------------------------------







#' read.nea.RData INPUT = Model Data (flows, inputs, outputs, storage)
#' formatted as for NEA.m, saved as CSV file S= |[F][z][X]| |[y][0][0]| OUPUT =
#' R Network data object for use with enaR
#' 
#' Borrett | July 15, 2013 --------------------------------------------------
#' read.nea.RData INPUT = Model Data (flows, inputs, outputs, storage)
#' formatted as for NEA.m, saved as CSV file S= |[F][z][X]| |[y][0][0]| OUPUT =
#' R Network data object for use with enaR
#' 
#' Borrett | July 15, 2013 --------------------------------------------------
#' Read NEA Formatted Network Model
#' 
#' This function reads in and creates a network object from a NEA formatted
#' data file (Fath and Borrett 2006).
#' 
#' @param file The name and path for the data file.
#' @param sep The separation character used to delimit data values.
#' @param warn LOGICAL: should pack warnings be reported?
#' @return Returns the network object.
#' @author Stuart R. Borrett
#' @seealso \code{\link{write.nea}}
#' @references Fath, B. D., Borrett, S. R. 2006. A Matlab function for Network
#' Environ Analysis.  Environ. Model. Softw. 21, 375-405.
#' @export read.nea
read.nea <- function(file="file name",sep=',',warn=TRUE){
  dat <- read.table(file,header=FALSE,sep=sep)  # assumes 
  n <- max(dim(dat)) - 2
  Flow <- t(dat[1:n,1:n])   # NEA.m stores flows col to row, so here we transpose
  z <- dat[1:n,(n+1)]  # inputs
  y <- dat[(n+1),1:n]  # outputs
  X <- dat[1:n,(n+2)]  # storage
  if (warn){
    model <- pack(flow=Flow,input=z,respiration=y,storage=X)  # create network data object
  }else{
    suppressWarnings(model <- pack(flow=Flow,input=z,respiration=y,storage=X))   # create network data object
  }
  return(model)
}
#' read.scor --- SCOR formatted file into R
#' in multiple formats
#' INPUT = file path
#' OUTPUT = network model in chosen format
#' S. Borrett and M. Lau | July 2011
#' ------------------------------------







#' read.scor --- SCOR formatted file into R in multiple formats INPUT = file
#' path OUTPUT = network model in chosen format S. Borrett and M. Lau | July
#' 2011 ------------------------------------ read.scor --- SCOR formatted file
#' into R in multiple formats INPUT = file path OUTPUT = network model in
#' chosen format S. Borrett and M. Lau | July 2011
#' ------------------------------------ Read SCOR Formatted Model
#' 
#' Read in network model data files that are in the SCOR format (REFERENCE).
#' 
#' The SCOR file must be formatted properly. In particular, the number of nodes
#' on the second line must have the first three characters dedicated to the
#' total number of nodes and the next three characters should contain the
#' number of living nodes. That is, the second line of the file should be
#' formatted as 'xxxyyy' where x and y are the characters for the total number
#' of nodes and the number of living nodes, respectively. Thus, if the total
#' number of nodes is 10 and the number of living nodes is 1, then the second
#' line should read, " 10 1."
#' 
#' @param file File path or plain text.
#' @param from.file States whether the file argument input should be treated as
#' a file path (TRUE) or plain text (FALSE).
#' @param warn Turn on (TRUE) or off (FALSE) warnings.
#' @return Returns the network model in one of several formats. The default
#' format is a network object used by the statnet package (type="network").
#' Three other options are the network environ analysis format (type="nea") as
#' defined by (Fath and Borrett 2006), a list format (type="list") and an edge
#' list (type="edge.list").
#' @author Matthew K. Lau Stuart R. Borrett
#' @references Ulanowicz, R.E. and J.J. Kay. 1991. A package for the analysis
#' of ecosystem flow networks. Environmental Software 6:131-142.
#' 
#' Fath, B. D., Borrett, S. R. 2006. A Matlab function for Network Environ
#' Analysis.  Environ. Model. Softw. 21, 375-405.
#' @export read.scor
read.scor <- function(file,from.file=TRUE,warn=FALSE){
  if (from.file){text <- readLines(file,warn=warn)}else{text <- file} # read in file
                                        #Partition the meta-data
  meta <- text[1]
                                        #Retrieve the number of node (n)
  n <- as.numeric(sub(' ','',substr(text[2],1,3)))
                                        #Determine the number of living nodes and create the living vector
  n.live <- as.numeric(sub(' ','',substr(text[2],4,6)))
  living <- c(rep(TRUE,n.live),rep(FALSE,(n-n.live)))
                                        #Retrieve vertex names
  vertex.names <- str_trim(text[3:(2+n)])
                                        # find negative ones (delimiters)
  br <- grep(pattern="( -1)", x=text)
  if (length(br) != 5){warning('Possible error in SCOR formatting')} # check expected number

  ## STORAGE
  B <- text[(3+n):(2+2*n)] # first cut at getting biomass data

  vertex.no <- as.numeric(sapply(B,function(x) (substr(x,1,3))))
  bm <- as.numeric(sapply(B,function(x) scifix(substr(x,5,nchar(x)))))
  storage <- data.frame("vertex"=vertex.no,"value"=bm) # final storage data

  ## INPUT

  if((br[2]-br[1])>1){
    inpt <- text[(br[1]+1):(br[2]-1)]
                                        #Condense into a function
    vertex.no <- as.numeric(sapply(inpt,function(x) substr(x,1,3)))
    z <- as.numeric(sapply(inpt,function(x) scifix(substr(x,5,nchar(x)))))
    inputs <- data.frame("vertex"=vertex.no,"value"=z) # final storage data

  } else {
    inputs <- NA
  }

  ## EXPORT

  if((br[3]-br[2])>1){
    export <- text[(br[2]+1):(br[3]-1)]
    vertex.no <- as.numeric(sapply(export,function(x) substr(x,1,3)))
    expt <- as.numeric(sapply(export,function(x) scifix(substr(x,5,nchar(x)))))
    exports <- data.frame("vertex"=vertex.no,"value"=expt) # final storage data

  } else {
    exports <- NA
  }


  ## RESPIRATION
  if((br[4]-br[3])>1){
    resp <- text[(br[3]+1):(br[4]-1)]
    vertex.no <- as.numeric(sapply(resp,function(x) substr(x,1,3)))
    resp <- as.numeric(sapply(resp,function(x) scifix(substr(x,5,nchar(x)))))
    respiration <- data.frame("vertex"=vertex.no,"value"=resp) # final storage data

  }else{
    respiration <- NA
  }


  ## FLOWS
                                        # assume there must be internal flows
  flows <- text[(br[4]+1):(br[5]-1)]
  strt <- as.numeric(sapply(flows,function(x) substr(x,1,3)))
  stp <- as.numeric(sapply(flows,function(x) substr(x,4,6)))
  value <- as.numeric(sapply(flows,function(x) scifix(substr(x,7,nchar(x)))))
  flows <- data.frame("tail"=strt,"head"=stp,"value"=value)
                                        #network data output type.
                                        # Note that the other output types depend on this sub-routine.

                                        #convert the flows to a matrix
  flow.mat <- array(0,dim=c(n,n))
  rownames(flow.mat) <- colnames(flow.mat) <- vertex.names
  for (i in seq(along=flows$tail)){
    flow.mat[flows$head[i],flows$tail[i]] <- flows$value[flows$tail==flows$tail[i]&flows$head==flows$head[i]]
  }
                                        #transpose flow matrix
  flow.mat <- t(flow.mat)
                                        #Vectorize the inputs
  input <- numeric(n)
  input[inputs$vertex] <- inputs$value
                                        #vectorize respiration and exports
  res <- numeric(n)
  exp <- numeric(n)

  if (any(is.na(exports)) == FALSE&any(is.na(respiration)) == FALSE){
    res[respiration$vertex] <- respiration$value
    exp[exports$vertex] <- exports$value
  }else if (any(is.na(exports))&any(is.na(respiration)) == FALSE){
    res[respiration$vertex] <- respiration$value
    exp <- rep(0,n)
  }else if (any(is.na(exports)) == FALSE&any(is.na(respiration))){
    exp[exports$vertex] <- exports$value
    res <- rep(0,n)
  }
  stor <- numeric(n)
  stor[storage$vertex] <- storage$value
                                        #produce an output vector given the values of export and respiration
  if (any(is.na(res)) == FALSE & any(is.na(exp))){
    output <- res
  }else if (any(is.na(res)) & any(is.na(exp)) == FALSE){
    output <- exp
  }else{output <- exp + res} #Outputs = respiration + exports for nea data type
                                        #introduce variables into the network format
  x <- pack(flow=flow.mat,input=input,export=exp,respiration=res,output=output,storage=stor,living=living)

  return(x)
}
#' read.wand --- WAND formatted file into R
#' INPUT = file path
#' OUTPUT = network object
#' S. Borrett | May 2012
#' ------------------------------------







#' read.wand --- WAND formatted file into R INPUT = file path OUTPUT = network
#' object S. Borrett | May 2012 ------------------------------------ read.wand
#' --- WAND formatted file into R INPUT = file path OUTPUT = network object S.
#' Borrett | May 2012 ------------------------------------ Read WAND Formatted
#' Model
#' 
#' Reads WAND formatted network models.
#' 
#' 
#' @param file File path to WAND formatted data file.
#' @return Returns a network object from a WAND formatted data file.
#' @note IMPORTANT: this function depends on the read.xlsx function from the
#' xlsx package, which requires that the entire path be specified from the root
#' directory (i.e. the absolute path).
#' @author Matthew K. Lau Stuart R. Borrett
#' @references Allesina, S., Bondavalli, C., 2004. WAND: an Ecological Network
#' Analysis user-friendly tool. Environmental Modelling and Software
#' 19(4):337-340.
#' @export read.wand
read.wand <- function(file='file name with path'){
                                        # file is the full excel file name
                                        # asssumes that first sheet is "Main" and second sheet is "Flows".
  x <- as.matrix(read.xls(file,sheet="Main"))
  d1 <- x[1:8,1] #model info
  n <- as.numeric(as.character(d1[3])) #Number of compartments
  dat.main <- x[8:(n+9),2:6] #isolate the stocks,imports,exports,respirations
  dat.main[is.na(dat.main)] <- 0 #zero NA
  vn <- dat.main[1:n,1] #vertex names
  dat.main <- apply(dat.main[,2:5],2,as.numeric)
                                        # get flows
  Flow <- read.xls(file,sheet="Flows")
  flow.mat <- as.matrix(Flow[1:(n),2:(n+1)])
  flow.mat[is.na(flow.mat)] <- 0
  flow.mat <- apply(flow.mat,2,as.numeric)
  rownames(flow.mat) <- colnames(flow.mat) <- vn
                                        #pack for export
  x <- list('flow'=flow.mat,
    'input'=dat.main[,2],
    'exports'=dat.main[,3],
    'respiration'=dat.main[,4],
    'storage'=dat.main[,1])
  ## --- Create Network Object From Data ---
  y <- network(x[[1]],directed=TRUE)
                                        # packing up the attributes into the network object (y)
  set.vertex.attribute(y,'input',x[[2]])
  set.vertex.attribute(y,'export',x[[3]])
  set.vertex.attribute(y,'respiration',x[[4]])
  set.vertex.attribute(y,'storage',x[[5]])
  set.vertex.attribute(y,'output',x[[3]]+x[[4]])
  y%v%'vertex.names' <- vn
  set.edge.attribute(y,'flow', flow.mat[flow.mat>0])
  return(y)
}
# relationalChange.r
# Which relationships change beteween the direct and integral utility analysis?
#
# Stuart R. Borrett
# February 25, 2016
# ------------------------







#' Relational change compared between two matrices.
#' 
#' Identifies the signs and pairwise relationsips of two matrices and compares
#' the difference between them.
#' 
#' 
#' @param x x is a square matrix of real numbers.  While this function is more
#' general, the initial intention was for this to be the direct utility matrix.
#' @param y y is a square matrix of real numbers.  While this function is more
#' general, the initial intention was for this to be the integral utility
#' matrix or the mixed trophic impacts matrix.
#' @return \item{Direct.Signs}{A sign matrix for matrix x.}
#' \item{Integral.Signs}{A sign matrix for matrix x.} \item{Direct.Relations}{A
#' matrix of the pairwise sign relationships for matrix x.}
#' \item{Integral.Relations}{A matrix of the pairwise signed relationships in
#' matrix y.} \item{Relations.Table}{A table that summarizes the relations.}
#' \item{Changed.Table}{A summary table of only the pariwise relationships that
#' changed between x and y.} \item{ns}{A vector of network statisitcs which
#' currently includes one whole-network statistic - a ratio of the
#' relationships changed between x and y.}
#' @note This function is called by enaUtility and enaMTI to summarize results.
#' @author Stuart R. Borrett
#' @seealso \code{\link{enaUtility}, \link{enaMTI}, \link{signs}}
#' @examples
#' 
#' 
#' 
#' data(oyster)
#' D <- enaUtility(oyster)$D
#' U <- enaUtility(oyster)$U
#' rc <- relationalChange(D, U)
#' 
#' 
#' ## To get a count of the number of differnt pairwise relationships in one of the
#' ## sign matrices, you can use the table function
#' 
#' count <- table(rc$Direct.Relations)
#' 
#' 
#' 
#' 
relationalChange <- function(x="Direct.U",y="Integral.U"){
    vnames <- rownames(x)
    S1 <- signs(x)    # find the signs of the relationships in the direct utility matrix
    S2 <- signs(y)    # find the signs of the relationships in the integral utility matrix
    S1$rs.tab$order <- 1:dim(S1$rs.tab)[1]  # add a column by which we can resort SF
    SF <- merge( S1$rs.tab, S2$rs.tab, by = c("From", "To"),stringsAsFactors=FALSE)  # merges the two relationship results
    names(SF) <- c("From","To","R1","R1.name","order", "R2", "R2.name")
    SF <- SF[,c(-4,-7)]  # remove relationship names (simplify)
    o <- order(SF$order)
    SF <- SF[o,!(names(SF) %in% c("order"))]  # reorder the merged data frame and drop order column

    # which pairs changed?
    d <- rep("-",dim(SF)[1])
    d[which(!(SF$R1 == SF$R2))] <- "*"  # find the differneces
    SF$changed <- d
    CR <- SF[which(SF$changed=="*"),]
    possible.change = length(d)
    no.changed = length(which(SF$changed=="*"))
    r.changed = round(no.changed/possible.change * 100, 2)

    ns <- c("r.change"=r.changed)  # percent of direct relationships that change when all utilities are considered

    return(list("Direct.Signs" = S1$s,
                "Direct.Relations" = S1$relations,
                "Integral.Signs" = S2$s,
                "Integral.Relations" = S2$relations,
                "Relations.Table" = SF,
                "Changed.Table" = CR,
                "ns"=ns))
}



#' scc --- find the strongly connected component
#' INPUT = an adjacency matrix
#' OUTPUT = list of membership and values
#' S. Borrett | July 2011
#' ------------------------------------







#' scc --- find the strongly connected component INPUT = an adjacency matrix
#' OUTPUT = list of membership and values S. Borrett | July 2011
#' ------------------------------------ scc --- find the strongly connected
#' component INPUT = an adjacency matrix OUTPUT = list of membership and values
#' S. Borrett | July 2011 ------------------------------------ Find the
#' Strongly Connected Component (SCC) in a Graph
#' 
#' This function finds the strongly connected components (SCCs) of an adjacency
#' matrix A and returns a number of derived network statistics.
#' 
#' 
#' @param A an n x n adjacency matrix.
#' @return \item{sp}{a list of structural properties including: the number of
#' SCCs ("no.scc"), the number of SCCs with more than 1 node ("no.scc.big"),
#' and the fraction of the network nodes participating in a large SCC ("pscc")}
#' \item{membership}{numeric vector giving the cluseter id to which each node
#' belongs (as in igraph:clusters)} \item{scc.id}{numeric vector of the numeric
#' identity in "membership" of SCCs with more than 1 node}
#' @note Input matrix is assumed to be oriented from columns to rows.
#' @author Matthew K. Lau Stuart R. Borrett
#' @seealso \link{enaStructure}
#' @references Allesina, S., Bodini, A., Bondavalli, C., 2005. Ecological
#' subsystems via graph theory: the role of strongly connected components.
#' Oikos 110, 164-176.
#' 
#' Berman, A., Plemmons, R.J., 1979. Nonnegative Matrices in the Mathematical
#' Sciences. Academic Press, New York.
#' 
#' Borrett, S.R., Fath, B.D., Patten, B.C. 2007. Functional integration of
#' ecological networks through pathway proliferation.  Journal of Theoretical
#' Biology 245, 98-111.
#' @examples
#' 
#' 
#' 
#' data(troModels)
#' A <- enaStructure(troModels[[6]])$A
#' scc(A)
#' 
#' 
#' 
#' @export scc
scc <- function(A="adjacency"){
                                        #Check for network class
  if (class(A) != 'matrix'){warning('A is not a matrix class object')}
  n <- dim(A)[1]
  c <- component.dist(A) # finds strong components in A (from sna package)
  no.scc <- length(c$csize)  # numer of scc
  j <- which(c$csize>1)  # finds scc > 1
  no.scc.big <- length(j)    # number of scc > 1
  pscc <- sum(c$csize[j])/n  # percent of nodes participating in a scc
  sp <- c("no.scc"=no.scc,"no.scc.big"=no.scc.big,"pscc"=pscc)
  y <- list("sp"=sp,"membership"=c$membership,"scc.id"=j-1)
  return(y)
}
#' scifix --- corrects missing e or E in 
#' scientific notation
#' INPUT = scalar either in or not in 
#' scientific notation
#' OUTPUT = corrected numeric value
#' M. Lau | July 2012
#' ------------------------------------







#' scifix --- corrects missing e or E in scientific notation INPUT = scalar
#' either in or not in scientific notation OUTPUT = corrected numeric value M.
#' Lau | July 2012 ------------------------------------ scifix --- corrects
#' missing e or E in scientific notation INPUT = scalar either in or not in
#' scientific notation OUTPUT = corrected numeric value M. Lau | July 2012
#' ------------------------------------ Standardizes Scientific Notation from
#' SCOR Formatted Files
#' 
#' This is a support function that corrects the scientific notation in SCOR
#' formatted data files.
#' 
#' 
#' @param x A numeric or character scalar.
#' @return Returns a numeric scalar in appropriate scientific notation.
#' @author Matthew K. Lau
#' @seealso \code{\link{read.scor}}
scifix <- function(x){
  x <- as.character(x)
                                        #e/E check
  e.check <- grepl('e',x)|grepl('E',x)
                                        #+/- check
  pm.check <- grepl('\\+',x)|grepl('\\-',x)
  if (pm.check&e.check==FALSE){
    if (grepl('\\+',x)){
      x <- sub('\\+','E+',x)
    }else{
      x <- sub('\\-','E-',x)
    }
  }else{}
  return(as.numeric(x))
}
#' set.orient --- globally reorients matrices
#' INPUT = matrix orientation (rc or cr)
#' OUTPUT = sets the expected orientation of matrices
#' 
#' M. Lau | Feb 2013
#' ---------------------------------------------------







#' set.orient --- globally reorients matrices INPUT = matrix orientation (rc or
#' cr) OUTPUT = sets the expected orientation of matrices
#' 
#' M. Lau | Feb 2013 ---------------------------------------------------
#' set.orient --- globally reorients matrices INPUT = matrix orientation (rc or
#' cr) OUTPUT = sets the expected orientation of matrices
#' 
#' M. Lau | Feb 2013 ---------------------------------------------------
#' Globally Set the Output Matrix Orientation
#' 
#' Changes the orientation of output matrices.
#' 
#' The enaR package as a whole, and the broader network analysis community,
#' assumes a row to column orientation; thus, the default orientation for the
#' package is row to column (DEFAULT = 'rc'). However, functions from the
#' Patten school were orignially developed to conduct calculations and produce
#' output in the column to row orientation. In order to facilitate the use of
#' these functions, we also provide the option for users to return output in
#' the orientation of the "school" (i.e. Patten results will be column to row
#' oriented) by setting the global orientation to "school" using this fuction.
#' 
#' @param x Orientation setting. If "rc" (DEFAULT), all matrix output will be
#' returned in row (=input) to column (=output) orientation, regardless of
#' school. If "school", then output matrices from functions from particular ENA
#' schools will be oriented as expected in that school (i.e. Patten =
#' column-row or Ulanowicz = row-column). Note, that all functions in the enaR
#' package expect input matrices to be oriented row-column.
#' @author Matthew K. Lau Stuart R. Borrett
#' @seealso \code{\link{get.orient}}
#' @examples
#' 
#' 
#' 
#' original.orientation = get.orient()
#' original.orientation
#' set.orient("school")
#' get.orient()
#' set.orient("rc")
#' get.orient()
#' set.orient(original.orientation)
#' 
#' 
#' 
#' @export set.orient
set.orient <- local({
  orientation <- 'rc'
  warn <- ''
  f <- function(x=c('rc','school')){
    if (any(x%in%c('rc','school','internal'))){
      orientation <<- x[1]
      if (x[1] == 'school'){
        warning('NOTE: output of functions from a particular analytical school will be returned in the standard orientation of that school.')
      }else if (x[1]=='internal'){
        orientation <<- 'school'
      }
    }else{
      warning('Unknown orientation.')
    }
  }
})
# signs.r
# February 25, 2016
# Stuart R. Borrett

# INPUT/OUPTUT - this function takes a numerical matrix and returns a matrix of its signs and a matrix of the relationships among the row/columns based on these signs
# -------------------







#' Signs and summary of input matrix
#' 
#' Identifies the signs and pairwise relationsips of a given matrix.  This
#' includes also returns a summary table that provides the ecological name of
#' each pairwise realtionship, and a summary of the counts.
#' 
#' 
#' @param x x is a square matrix of real numbers.  While this function is more
#' general, the initail intention was for this to be a utility matrix or the
#' mixed trophic impacts matrix.
#' @return \item{sign}{A sign matrix for matrix x.} \item{relations}{A matrix
#' of the pairwise signed relationships in x.} \item{rs.tab}{Table summarizing
#' the pairwise relationships and identifying their ecological label.}
#' \item{relationship.counts}{A count of the different kinds of pairwise
#' relationships found in matrix x.}
#' @note This function is called by relationalChange, and was created to
#' generate more informative output from enaUtility and enaMTI.
#' @author Stuart R. Borrett
#' @seealso \code{\link{relationalChange}}
#' @examples
#' 
#' 
#' 
#' data(oyster)
#' U <- enaUtility(oyster)$U
#' s <- signs(U)
#' 
#' 
#' 
signs <- function(x="matrix"){
    vnames <- rownames(x)   # get row names
    d <- dim(x)             # find matrix dimensions
    s <- matrix(rep(0,prod(d)), ncol = d[2])  # initialize sign matrix
    r <- s                                    # initialize relationship matrix
    r.sparse <- list()                        # initialize r.sparse
    cnt1 <- 1                                 # iniitlize counter

    # dictionary of ecological relationship types
    erd <- matrix(c( "-", "+", "(-,+)",  "predation",
                    "-", "0", "(-,0)", "amensalism",
                    "-", "-", "(-,-)",  "competition",
                    "0", "+", "(0,+)", "anabolism",
                    "0", "0", "(0,0)", "neutralism",
                    "0", "-", "(0,-)", "catabolism",
                    "+", "+", "(+,+)", "mutualism",
                    "+", "0", "(+,0)", "commensalism",
                    "+", "-", "(+,-)", "altruism"),
                  ncol  = 4, byrow = TRUE)
    colnames(erd) <- c("source","sink","pair","name")

    # FIND matrix element sign {+, 0, -}
    positive <- which(x > 0, arr.ind = TRUE)  # which are positive
    negative <- which(x < 0, arr.ind = TRUE)  # which are negative

    # build the sign matrix
    s[positive] <- "+"
    s[negative] <- "-"
    rownames(s) <- colnames(s) <- vnames

    # combine signs to determine pairwise relationships
    for(i in 1:d[1]){
        for(j in 1:d[2]){
            if(j>=i){
                r[i,j] <- paste("(",s[j,i],",",s[i,j],")",sep="")  # pariwise realtionships - matrix
                tmp.relation <- erd[which(erd[,3] == r[i,j]),4]
                r.sparse[[cnt1]] <- c(vnames[i],vnames[j],r[i,j],tmp.relation)  # pariwise realtionships - table
                cnt1 <- cnt1 + 1  # increment counter
            }
        }
    }
    colnames(r) <- vnames; rownames(r) <- vnames
    r.sparse <- do.call(rbind,r.sparse)  # combine results for r.sparse
    colnames(r.sparse) <- c("From","To","Relationship","R.name")  # rename columns

    relationship.counts <- table(r.sparse[,3])

    return(list("sign"=s,
                "relations"=r,
                "rs.tab"=as.data.frame(r.sparse,stringsAsFactors=FALSE),
                'relationship.counts' = relationship.counts)
           )
}
#' ssCheck --- checks if the given network
#' is out of balance by a given tolerance
#' threshold
#' INPUT = network object
#' OUTPUT = logical indicating violation of
#' tolerance
#' NOTE: used in the balancing process
#' M. Lau | July 2011
#' ------------------------------------







#' ssCheck --- checks if the given network is out of balance by a given
#' tolerance threshold INPUT = network object OUTPUT = logical indicating
#' violation of tolerance NOTE: used in the balancing process M. Lau | July
#' 2011 ------------------------------------ ssCheck --- checks if the given
#' network is out of balance by a given tolerance threshold INPUT = network
#' object OUTPUT = logical indicating violation of tolerance NOTE: used in the
#' balancing process M. Lau | July 2011 ------------------------------------
#' Checks the Balance of Inputs and Outputs from a Network
#' 
#' This function supports the balancing process by checking if the inputs and
#' outputs of a given network model are within acceptable limits.
#' 
#' 
#' @param x A network object.
#' @param tol The threshold for balance in percent difference between input and
#' outputs.
#' @param more LOGICAL: should more detailed results be returned?
#' @param zero.na LOGICAL: should NA values be changed to zeros?
#' @return Returns a logical value stating if the model is within acceptable
#' limits of balance (TRUE) or if it is not (FALSE).
#' @author Matthew K. Lau Stuart R. Borrett
#' @seealso \code{\link{balance}}
#' @references Fath, B.D. and S.R. Borrett. 2006. A MATLAB function for network
#' environ analysis. Environmental Modelling & Software 21:375-405.
#' @examples
#' 
#' 
#' 
#' data(troModels)
#' ssCheck(troModels[[2]])
#' ssCheck(troModels[[6]])
#' 
#' 
#' 
#' @export ssCheck
ssCheck <- function(x,tol=5,more=FALSE,zero.na=TRUE){
                                        #Check for network class object
  if (class(x) != 'network'){warning('x is not a network class object')}
  T. <- as.extended(x) #convert to extended format
  if (zero.na){T.[is.na(T.)] <- 0}
  n <- network.size(x) #get the number of nodes
  Tin <- apply(T.[,1:n],2,sum) #in throughflow
  Tout <- apply(T.[1:n,],1,sum) #out throughflow
  d <- abs(Tin - Tout) # SSerror difference
  pe <- (d / Tout)*100 # SSerror as percent of total throughflow
                                        #
  if(more==FALSE){
    return(all(pe < tol)) #returns a logical indicating that all node differences are less than tolerance (==TRUE)
  }else{
    return(list("ss"=all(pe < tol),"Tin"=Tin,"Tout"=Tout,"perror"=pe))
  }
}
#' structure.statistics --- calculates structural statistics
#' INPUT = an adjacency matrix
#' OUTPUT = list of structural statistics
#' S. Borrett | July 2011
#' ------------------------------------







#' structure.statistics --- calculates structural statistics INPUT = an
#' adjacency matrix OUTPUT = list of structural statistics S. Borrett | July
#' 2011 ------------------------------------ structure.statistics ---
#' calculates structural statistics INPUT = an adjacency matrix OUTPUT = list
#' of structural statistics S. Borrett | July 2011
#' ------------------------------------ Structural Statistics of an Ecological
#' Network
#' 
#' This function returns several network statistics that describe a network.
#' 
#' 
#' @param A An adjacency matrix.
#' @return \item{n}{Number of nodes in A.} \item{L}{Number of direct
#' connections in A.} \item{C}{Connectivity of A.} \item{LD}{Link density.}
#' \item{lam1A}{First dominant eigenvalue of A.} \item{mlam1A}{Multiplicity of
#' the dominant eigenvalue.} \item{lam2A}{Magnitude of the second largest
#' eigenvalue.} \item{rho}{Damping ratio (see Caswell 2001).} \item{R}{Distance
#' of lam1A from the bulk of the eigen spectrum.} \item{d}{Difference between
#' the dominant eigenvalue and the link density.} \item{no.scc}{Number of
#' strongly connected components.} \item{no.scc.big}{Number of strongly
#' connected components greater than 1.} \item{pscc}{Precent of nodes
#' participating in a strongly connected component.}
#' @author Matthew K. Lau Stuart R. Borrett
#' @seealso \code{\link{enaStructure}},\code{\link{scc}}
#' @references Fath, B. D., Borrett, S. R. 2006. A Matlab function for Network
#' Environ Analysis.  Environ. Model. Softw. 21, 375-405.
structure.statistics <- function(A='adjacency matrix'){
  if (class(A) != 'matrix'){warning('A is not a matrix class object')}
                                        # 
  n <- dim(A)[1] #number of nodes in A
  L <- sum(A)    #length(A[A!=0]) #number of direct connections in A
  C <- L/n^2 #connectivity of A
  LD <- L/n #link density (equivalent to n*C)
  e <- eigen(A)$values
  aer <- round(abs(e),digits = 7)      # round eigenvalue magnitudes (remove numerical error)
  mlam1A <- length(which(aer == aer[1]))  # find multiplicity of dominant eigenvalue
  ppr <- sum(mExp(A,200))/sum(mExp(A,199)) # pathway proliferation rate (Borrett & Patten 2003)
  lam1A <- abs(e[1])                   # dominant eigenvalue of A. Also termed spectral radius. This is 1) a measure of connectivity, 2) approximately equal to LD and 3) the rate of pathway proliferation
  d <- abs(lam1A-LD)                   # difference between dominant eigenvalue and link density
  
  if ((n-mlam1A)>0){
    lam2A <- abs(e[(1+mlam1A)]) #magnitude of second largest eigenvalue
    rho <- lam1A/abs(lam2A)  #damping ratio, an indicator of how quickly a^(m)/a^(m-1) foes to lam1[A] (Caswell 2001, p. 95)
    R <- abs(e[n])-abs(e[n-1])/(abs(e[n-1])-abs(e[1])) #distance of lam1[A] from the bulk of the eigen spectrum (Farkas et al. 2001)
  }  else {
    lam2A <- NA;rho <- NA;R <- NA #flag to indicate these do not exist
  }	
  sp1 <- as.vector(scc(A)$sp)
  no.scc <- sp1[1]
  no.scc.big <- sp1[2]
  pscc <- sp1[3]
  sp <- cbind(n,L,C,LD,ppr,lam1A,mlam1A,rho,R,d,no.scc,no.scc.big,pscc)  # list of structural statistics of interest
  return(sp)
}
#' unpack --- extracts network object into
#' a list
#' INPUT = network object
#' OUTPUT = list of network model components
#' S. Borrett and M. Lau | July 2011
#' ------------------------------------







#' unpack --- extracts network object into a list INPUT = network object OUTPUT
#' = list of network model components S. Borrett and M. Lau | July 2011
#' ------------------------------------ unpack --- extracts network object into
#' a list INPUT = network object OUTPUT = list of network model components S.
#' Borrett and M. Lau | July 2011 ------------------------------------
#' "Unpacks" the Network Object into Separate Objects
#' 
#' Separates the components of a network object into separate components within
#' a list. This includes inputs, exports, respirations, outputs (exports +
#' respirations), storage, and internal flows.
#' 
#' 
#' @param x A network object.  This includes all weighted flows into and out of
#' each node.
#' @return \item{F}{matrix of flows from each node to each node oreinted row to
#' column.} \item{z}{Node boundary inputs.} \item{r}{Node boundary loss from
#' respiration.} \item{e}{Node boundary loss due to exportation} \item{y}{Node
#' boundary loss; summation of r and e} \item{X}{Node storage or biomass}
#' \item{living}{Logical vector indicating whether each node is living or not}
#' @note Flows are oriented from row to column.
#' @author Matthew K. Lau Stuart R. Borrett
#' @seealso \code{\link{pack},\link{read.scor}}
#' @examples
#' 
#' 
#' 
#' data(troModels)
#' unpack(troModels[[6]])
#' 
#' 
#' 
#' @export unpack
unpack <- function(x='network object'){
  flow <- as.matrix(x, attrname = 'flow')
  input <- x%v%'input'
  respiration <- x%v%'respiration'
  respiration[is.na(respiration)] <- 0
  export <- x%v%'export'
  export[is.na(export)] <- 0
  output <- x%v%'output'   #respiration + export
  storage <- x%v%'storage'
  living <- x%v%'living'
  return(list("F"=flow,"z"=input,"r"=respiration,"e"=export,"y"=output,"X"=storage,'living'=living))
}
###Function to output a model for EcoNet
###http://eco.engr.uga.edu/DOC/econet2.html
###MKLau 17nov2014







#' Write enaR models to an EcoNet formatted file.
#' 
#' Creates an EcoNet model from an enaR network object that can be used with
#' the online interface for EcoNet.
#' 
#' 
#' @param x Network object.
#' @param file The file name or path. If a simple file name is given, this
#' function uses the current working directory by default.
#' @param mn The model name that EcoNet will use. The DEFAULT is 'ena_model'.
#' @param zero.flows LOGICAL: should zero flow values be written?
#' @return An EcoNet formatted text file is created from the model, which can
#' be input at http://eco.engr.uga.edu.
#' @author Matthew K. Lau
#' @references About EcoNet (http://eco.engr.uga.edu/DOC/econet1.html)
#' 
#' Kazanci, C. 2009. Handbook of Ecological Modelling and Informatics, by WIT
#' Press.
#' @export write.EcoNet
write.EcoNet <- function(x='model',file='file path',mn='ena_model',zero.flows=FALSE){
    x <- unpack(x)
###node names
    nn <- rownames(x$F)
    nn <- strsplit(nn,split='')
    nn <- lapply(nn,function(x) x[x%in%letters|x%in%LETTERS|x%in%(1:9)])
    nn <- unlist(lapply(nn,paste,collapse=''))
    rownames(x$F) <- colnames(x$F) <- nn
###Write model name
    write.table(paste('###',mn),file=file,col.names=FALSE,row.names=FALSE,quote=FALSE)
###initial conditions
    write.table(paste('###','initial conditions'),file=file,col.names=FALSE,row.names=FALSE,quote=FALSE,append=TRUE)
    write.table(paste(rownames(x$F),x$X,sep='='),file=file,col.names=FALSE,row.names=FALSE,quote=FALSE,append=TRUE)
###write inputs
    write.table(paste('###','inputs'),file=file,col.names=FALSE,row.names=FALSE,quote=FALSE,append=TRUE)
    write.table(paste('*',' -> ',rownames(x$F),' c=',x$z,sep=''),file=file,col.names=FALSE,row.names=FALSE,quote=FALSE,append=TRUE)
###write flows
    write.table(paste('###','flows'),file=file,col.names=FALSE,row.names=FALSE,quote=FALSE,append=TRUE)
    for (i in 1:nrow(x$F)){
        for (j in 1:nrow(x$F)){
            if (x$F[i,j] == 0){}else{
                write.table(paste(rownames(x$F)[i],' -> ',colnames(x$F)[j],' c=',x$F[i,j],sep=''),
                            file=file,col.names=FALSE,row.names=FALSE,
                            quote=FALSE,append=TRUE)
            }
        }
    }
###write outputs
    write.table(paste('###','outputs'),file=file,col.names=FALSE,row.names=FALSE,quote=FALSE,append=TRUE)
    write.table(paste(rownames(x$F),' -> ','*',' c=',(x$r+x$e),sep=''),file=file,col.names=FALSE,row.names=FALSE,quote=FALSE,append=TRUE)
}
#'' write.nea.R
#' INPUT = enaR network data object
#' Ouput = CSV formatted file with data arranged as expected input for NEA.m
#'
#' Borrett | July 15, 2013
#' ----------------------------------------







#' ' write.nea.R INPUT = enaR network data object Ouput = CSV formatted file
#' with data arranged as expected input for NEA.m
#' 
#' Borrett | July 15, 2013 ---------------------------------------- '
#' write.nea.R INPUT = enaR network data object Ouput = CSV formatted file with
#' data arranged as expected input for NEA.m
#' 
#' Borrett | July 15, 2013 ---------------------------------------- Write a
#' Network Object to File Using the NEA Data Format
#' 
#' This function writes a network object to a NEA formatted data file (Fath and
#' Borrett 2006).
#' 
#' @param x Network object.
#' @param file.name The file name or path. If a simple file name is given, this
#' function uses the current working directory by default.
#' @param sep The separation character used to delimit data values.
#' @return Writes a network object to a NEA formatted file and returns the
#' output composite matrix.
#' @author Stuart R. Borrett
#' @seealso \code{\link{read.nea}}
#' @references Fath, B. D., Borrett, S. R. 2006. A Matlab function for Network
#' Environ Analysis.  Environ. Model. Softw. 21, 375-405.
#' @export write.nea
write.nea <- function(x, file.name,sep=','){
                                        # Check for network class
  if (class(x) != 'network'){warning('x is not a network class object')}
  U <- unpack(x)  # unpack data
  n <- length(U$z)
  S <- matrix(NA,nrow=(n+1),ncol=(n+2))
  S[1:n,1:n] = t(U$F)
  S[1:n,(n+1)]= U$z
  S[(n+1),1:n] = U$y
  S[1:n,(n+2)]= U$X
  S[(n+1),(n+1):(n+2)] = 0
                                        # write file
  write.table(S,file=file.name,row.names=FALSE,col.names=FALSE,sep=sep) 
                                        # return composite system matrix to workspace 
  return(S)
}
