# Example enaR Code
# mee13 manuscirpt
# October 11, 2013
# ---------------------
rm(list=ls())
library(enaR)   # load enaR package

# -- ENTER MODEL DATA -- from Dame and Patten (1981)
# node names
names <- c("Filter Feeders","Microbiota","Meiofauna",
                      "Deposit Feeders","Predators","Deposited Detritus")

# Internal Flows of model, as matrix (oriented row to column)
F <- matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 8.1721, 0, 1.2060, 0, 0, 0, 7.2745,
               0, 1.2060, 0.6609, 0, 0, 0.6431, 0.5135, 0, 0,
               0.1721, 0, 0, 15.7910, 0, 4.2403, 1.9076, 0.3262, 0), ncol=6)
rownames(F) <- names  # add node names to rows
colnames(F) <- names  # add node names to cols

# boundary flows
inputs <- c(41.47,0, 0, 0, 0, 0)
outputs <- c(25.1650, 5.76, 3.5794, 0.4303, 0.3594, 6.1759)

# Living
Living <- c(TRUE,TRUE,TRUE,TRUE,TRUE,FALSE)

# pack the model data into the R network data object
m <- pack(flow=F,input=inputs, respiration=outputs,output=outputs, living=Living)
ssCheck(m)

### THIS IS NO LONGER NEEDED
# m <- read.scor("oyster.dat")  # read model data from SCOR formatted file
# m <- balance(m)               # balance model using AVG2 algorithm (Allesina and Bondavali)
# u <- unpack(m)                # unpack model data to illustrate components
# attributes(u)

F <- enaFlow(m)               # perform ENA flow analysis
attributes(F)                  # show analysis objects created
F$ns                          # show flow analysis network statistics
F$T

