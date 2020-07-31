`%+%` <- function(x, y)  mapply(sum, x, y, MoreArgs = list(na.rm = TRUE)) #skip NAs when adding vectors

sbml2dnf <- function(SBML){
  "Converts the regulatory functions in an SBML objects attribute 'interactions' into their disjunctive normal forms. 
  Used to creative trinary adjacency matrices with trinaryadjmat(SBML). 
  For binary adjacency matrices, use only the SBML object in its unmodified form."
  saveNetwork(SBML, file = "tmpnet.txt", generateDNFs = T)
  SBML_DNF <- loadNetwork("tmpnet.txt", symbolic = T)
  file.remove("tmpnet.txt")
  return(SBML_DNF)
}

graph_to_amatrix <- function(g) {
  "Converts an igraph graph object to an adjacency matrix"
  a <- get.adjacency(g)
  mat <- matrix(as.numeric(a), ncol=ncol(a))
  for (d in 1:dim(mat)[2]){if (mat[d,d] >= 1){mat[d,d] <- 1}}
  return(mat)
}

conv2adjmat <- function(sbmlnet, inputcorrected = FALSE){
  "Converts an SBML object to adjacency matrix. 
  If inputcorrected = T, all input edges are set to zero. 
  A vertex v is said to be an input if it is only regulated by itself, 
  meaning the sum of column v in the adjacency matrix is one, with the only non-zero entry 
  being at position [v,v]."
  adjmat <- sapply(sbmlnet$interactions, function(gene) {v <- rep(0,length(sbmlnet$genes)); 
  v[gene$input] <- 1; return(v)})
  if (inputcorrected == TRUE){
    for (d in 1:dim(adjmat)[1]){
        if (adjmat[d,d] >= 1){adjmat[d,d] <- 1}
      }
    for (col in 1:dim(adjmat)[2]){
      if (sum(adjmat[,col]) == 1 & adjmat[col,col] == 1){
        adjmat[col,col] <-0
      }
    }
  }
  return(adjmat)
}

totalgenes <- function(adjmat){
  #Returns total number of vertices (i.e. genes)
  return(dim(adjmat)[1])
}

generateuniquestates <- function(N=2^l, l, fixGenes, returnList = T){
  "Returns a matrix with N rows, each row being a random unique Boolean sequence of length l.
  fixGenes is an optional vector of states that should have a certain value. Needs to be of length l.
  If fixGenes[g] = 1 or 0 -> remove all states where grid[,g] is 0 or 1. If fixGenes[g] has any other value, ignore it."
  if (N > 2^l){return(paste("N needs to be smaller than 2^l=", 2^l, sep=""))}
  if (missing(fixGenes) ==F & missing(N) == T){N <- 2^(l-length(which(fixGenes %in% c(1,0) ==T)))}
  #return matrix with N random unique binary vectors of length l
  states <- rep(list(0:1), l)
  grid <- expand.grid(states)
  colnames(grid) <- seq(1,l)
  if (N==2 & l==1){
    if (returnList == F){
      return(grid)
    } else{
      return(lapply(seq_len(nrow(grid)), function(i) grid[i,]))
    }
  }
  
  if (missing(fixGenes)){
    randomstates <- sample(1:dim(grid)[1], N, replace=F)
    grid <- grid[randomstates,]
    rownames(grid) <- seq(1, dim(grid)[1])
    if (returnList == F){
      return(grid)
    } else {
      return(lapply(seq_len(nrow(grid)), function(i) grid[i,]))
    }
  }
  else {
    nrfixed <- length(which(fixGenes %in% c(1,0) ==T))
    maxpossible <- 2^(l-nrfixed)
    if (length(fixGenes) != l){return(paste("Vector fixGenes must have length l=", l, sep=""))}
    for (g in 1:l){
      if (fixGenes[g] == 1){grid <- grid[which(grid[,g]==1),]}
      else if (fixGenes[g] == 0){grid <- grid[which(grid[,g]==0),]}
    }
    if (dim(grid)[1] < N){return(paste("Error: Not enough states left. Must have N<=", maxpossible, ".", sep=""))}else{
      randomstates <- sample(1:dim(grid)[1], N, replace=F)
      grid <- grid[randomstates,]
      rownames(grid) <- seq(1, dim(grid)[1])
    }
    if (returnList == F){
      return(grid)
    } else {
      return(lapply(seq_len(nrow(grid)), function(i) grid[i,]))
    }
  }
}

VertexBetweenness <- function(adjmat, normalized=T){
  " Returns the vertex betweenness for all n vertices in the network. 
  If normalized=T, each value is normalised by the factor 1/(n*n-3*n+2)."
  graph <- graph_from_adjacency_matrix(adjmat)
  #returns double with betweenness centrality for every vertex
  return(betweenness(graph, normalized = normalized))
}

flattenAttractor <- function(attractor, nr, equidistant=F){
  "Takes an attractor object from the BoolNet function getAttractors() and the index of a specific attractor 
  in this object. Returns a trinary vector whose length corresponds to the total number of genes in 
  the network where each value indicates a pattern of gene expression. 
  '1' and '0' indicate the gene being always turned on or always turned off in the attractor states respectively, 
  while '2' indicates any kind of oscillation. If \texttt{equidistant=T}, 
  the value of oscillating genes is set to '0.5' instead."
  totalattrs <- length(attractor[[2]])
  if (nr > totalattrs){return(paste("Error: There are only", totalattrs, "attractors in total."))}
  length <- length(attractor$stateInfo$genes)
  #Decimal zu Bin?r
  attr <- apply(attractor[[2]][[nr]]$involvedStates, MARGIN = 2, BoolNet:::dec2bin, len = length)
  #fixed oder oszillierend?
  fixedT <- apply(attr, MARGIN = 1, function(x) all(x == T))
  fixedF <- apply(attr, MARGIN = 1, function(x) all(x == F))
  
  attractorState <- rep(2, length)
  
  attractorState[which(fixedT)] <- 1
  attractorState[which(fixedF)] <- 0
  if (equidistant==F){return(attractorState)}
  else{
    return(replace(attractorState, attractorState==2, 0.5))
  }
  #modify to return correct result for attractor of length 1 (necessarily only 1 or 0)
}
probactivation <- function(attractor, nr, gene){
  totalattrs <- length(attractor[[2]])
  if (nr > totalattrs){return(paste("Error: There are only", totalattrs, "attractors in total."))}
  length <- length(attractor$stateInfo$genes)
  attr <- apply(attractor[[2]][[nr]]$involvedStates, MARGIN = 2, BoolNet:::dec2bin, len = length)
  genets <- attr[gene,] #time series of expression of GENE in given attractor NR
  probofexpr <- sum(genets)/length(genets)
  return(probofexpr)
}
probAttractor <- function(attractor, nr){
  "Takes an attractor object obtained using the BoolNet function getAttractors() and an attractor index. 
  Returns the probability for any given gene to be active during an attractor cycle."
  total <- length(attractor$attractors)
  dim <- length(attractor$stateInfo$genes)
  res <- rep(NA, dim)
  for (g in 1:dim){
    res[g] <- probactivation(attractor, nr, g)
  }
  return(res)
}

#Dynamic measure - Hamming distance
meanofminattrdists <- function(SBML, modification, modgene, method = "trinary"){
  if (method %in% c("trinary", "prob") == F){return("Not a valid method, choose trinary or prob.")}
  if (modification %in% c("overexpression", "knockout") == F){return("Not a valid modification, choose overexpression or knockout.")}
  attrs <- getAttractors(SBML, method = "sat.exhaustive")
  modSBML <- SBML
  if (modification == "overexpression"){modSBML$fixed[modgene] <- 1}
  else if (modification == "knockout"){modSBML$fixed[modgene] <- 0}
  modattrs <- getAttractors(modSBML, method = "sat.exhaustive")
  res <- rep(NA, length(modattrs$attractors))
  for (ma in 1:length(modattrs$attractors)){
    if (method == "trinary"){
      flattenedmod <- flattenAttractor(modattrs, ma, equidistant = T)
    } else {flattenedmod <- probAttractor(modattrs, ma)}
    flattenedmod <- flattenedmod[-modgene]
    dists <- rep(NA, length(attrs$attractors))
    for (a in 1:length(attrs$attractors)){
      if (method == "trinary"){
        flattenedorig <- flattenAttractor(attrs, a, equidistant = T)
      } else {flattenedorig <- probAttractor(attrs, a)}
      flattenedorig <- flattenedorig[-modgene]
      dist_ma_from_a <- sum(abs(flattenedorig - flattenedmod))
      dists[a] <- dist_ma_from_a
    }
    res[ma] <- min(dists)
  }
  return(mean(res))
}
meanofminattrdistsbygene <- function(SBML, modification, method = "trinary"){
  "Modifies every gene either by overexpression or knockout and returns the mean of the minimal distances 
  between a fixed modified and all unmodified attractors as measured by their overlaps, given either the 
  method 'trinary' or 'prob'."
  if (modification %in% c("overexpression", "knockout") == F){return("Not a valid modification, choose overexpression or knockout.")}
  res <- rep(NA, length(SBML$genes))
  for (g in 1:length(SBML$genes)){
    res[g] <- meanofminattrdists(SBML, modification, g, method)
  }
  return(res)
}
meanofminattrdistsmerged <- function(SBML, method = "trinary"){
  "Returns a vector choosing the worst-case output from either OE or KO for any gene, 
  i.e. the most impactful perturbation."
  res <- OE <- KO <- rep(NA, length(SBML$genes))
  for (g in 1:length(SBML$genes)){
    OE[g] <- meanofminattrdists(SBML, modification="overexpression", g, method)
    KO[g] <- meanofminattrdists(SBML, modification="knockout", g, method)
    res[g] <- max(OE[g], KO[g])
  }
  names(res) <- SBML$genes
  return(res)
}

#Dynamic measure - Attractor loss
attrcoverageaftermod <- function(SBML, modification, modgene, method = "trinary"){
  if (method %in% c("trinary", "prob") == F){return("Not a valid method, choose trinary or prob.")}
  if (modification %in% c("overexpression", "knockout") == F){return("Not a valid modification, choose overexpression or knockout.")}
  attrs <- getAttractors(SBML, method = "sat.exhaustive")
  modSBML <- SBML
  if (modification == "overexpression"){modSBML$fixed[modgene] <- 1}
  else if (modification == "knockout"){modSBML$fixed[modgene] <- 0}
  modattrs <- getAttractors(modSBML, method = "sat.exhaustive")
  #Check how many orig. attractors have a matching modattr. with distance = 0 (phenotype preserved, excluding modgene)
  losscounter <- 0 #count phenotypes that have no match in modSBML, finally divide by nr. of orig.attrs
  for (a in 1:length(attrs$attractors)){
    dists <- rep(NA, length(modattrs$attractors))
    if (method == "trinary"){
      flattenedorig <- flattenAttractor(attrs, a, equidistant = T)
    } else {flattenedorig <- probAttractor(attrs, a)}
    flattenedorig <- flattenedorig[-modgene]
    for(ma in 1:length(modattrs$attractors)){
      if (method == "trinary"){
        flattenedmod <- flattenAttractor(modattrs, ma, equidistant = T)
      } else {flattenedmod <- probAttractor(modattrs, ma)}
      flattenedmod <- flattenedmod[-modgene]
      dists[ma] <- sum(abs(flattenedorig-flattenedmod))
    }
    if (min(dists) > 0){losscounter <- losscounter + 1}
  }
  return(1 - (losscounter/length(attrs$attractors)))
}
attrcoverageaftermodbygene <- function(SBML, modification, method = "trinary"){
  "Modifies every gene either by overexpression or knockout and returns the percentage of original attractors 
  that are still existent in the modified network, that is those with a minimal distance of 0 to a modified 
  attractor, as measured by their overlaps, given either the \texttt{method} 'trinary' or 'prob'."
  if (modification %in% c("overexpression", "knockout") == F){return("Not a valid modification, choose overexpression or knockout.")}
  res <- rep(NA, length(SBML$genes))
  for (g in 1:length(SBML$genes)){
    res[g] <- attrcoverageaftermod(SBML, modification, g, method)
  }
  names(res) <- SBML$genes
  return(res)
}
attrcoveragemerged <- function(SBML, merging="min", method="trinary"){
  if (merging %in% c("min", "mean") == F){return("Error: merging must be 'min' or 'mean'.")}
  if (method %in% c("trinary", "prob") == F){return("Error: method must be 'trinary' or 'prob'.")}
  #merging=c("mean", min") - mean: average OE + KO coverage; min - worst case, pick lower coverage
  #method=c("trinary", "prob")
  res <- OE <- KO <- rep(NA, length(SBML$genes))
  for (g in 1:length(SBML$genes)){
    OE[g] <- attrcoverageaftermod(SBML, modification="overexpression", g, method=method)
    KO[g] <- attrcoverageaftermod(SBML, modification="knockout", g, method=method)
    if (merging=="mean"){res[g] <- (OE[g]+KO[g])/2}
    else if (merging=="min"){res[g] <- min(OE[g], KO[g])}else{return("Not a valid merging scheme.")}
  }
  names(res) <- SBML$genes
  return(res)
}

#Dynamic measure - Attractor gain
newattrgain <- function(SBML, modification, modgene, method="trinary"){
  if (method %in% c("trinary", "prob") == F){return("Not a valid method, choose trinary or prob.")}
  if (modification %in% c("overexpression", "knockout") == F){return("Not a valid modificiation!")}
  #mod every gene, count occurrences of NEW attractors, i.e. distance to ALL orig attrs is >0
  adjmat <- conv2adjmat(SBML, inputcorrected = F)
  size <- dim(adjmat)[1]
  attrs <- getAttractors(SBML, method = "sat.exhaustive")
  nrattrs <- length(attrs$attractors)
  modSBML <- SBML
  if (modification=="overexpression"){modSBML$fixed[modgene] <- 1} else if (modification=="knockout"){modSBML$fixed[modgene] <- 0}
  modattrs <- getAttractors(modSBML, method = "sat.exhaustive")
  nrmodattrs <- length(modattrs$attractors)
  res <- nrmodattrs #start by assuming all modattrs are new, -1 if a modattr has an orig. equivalent with distance =0
  #Check how many orig. attractors have a matching modattr. with distance = 0 (phenotype preserved, excluding modgene)
  for(ma in 1:nrmodattrs){
    dists <- rep(NA, length(attrs$attractors)) #distances of given ma to all available a
    if (method == "trinary"){
      flattenedmod <- flattenAttractor(modattrs, ma, equidistant = T)
    } else {flattenedmod <- probAttractor(modattrs, ma)}
    flattenedmod <- flattenedmod[-modgene]
    for (a in 1:nrattrs){
      if (method == "trinary"){
        flattenedorig <- flattenAttractor(attrs, a, equidistant = T)
      } else {flattenedorig <- probAttractor(attrs, a)}
      flattenedorig <- flattenedorig[-modgene]
      dists[a] <- sum(abs(flattenedorig-flattenedmod))
    }
    if (min(dists) == 0){res <- res - 1} #there exists an orig attr with dist=0 to ma -> ma is not a newly gained attr
  }
  return(res)
}
newattrgainbygene <- function(SBML, modification, method="trinary"){
  adjmat <- conv2adjmat(SBML, inputcorrected = F)
  size <- dim(adjmat)[1]
  res <- rep(NA, size)
  for (g in 1:size){
    res[g] <- newattrgain(SBML=SBML, modification=modification, modgene = g, method=method)
  }
  names(res) <- SBML$genes
  return(res)
}
newattrgainbygenemerged <- function(SBML, method="trinary"){
  adjmat <- conv2adjmat(SBML, inputcorrected = F)
  size <- dim(adjmat)[1]
  res <- rep(NA, size)
  OEimpact <- newattrgainbygene(SBML, modification="overexpression", method=method)
  KOimpact <- newattrgainbygene(SBML, modification="knockout", method=method)
  for (g in 1:size){
    res[g] <- max(OEimpact[g], KOimpact[g])
  }
  names(res) <- colnames(adjmat)
  return(res)
}

#Static measure - Determinative Power
binaryShannonEntropy <- function(p){
  "Get binary Shannon entropy h(p), probability p must be [0,1]."
  if (p < 0 | p > 1){return("Error: Probability p must be in range [0,1].")}
  if (p == 0 | p == 1){return(0)}
  return(-p*log2(p)-(1-p)*log2(1-p))
}
suppBoolFunc <- function(SBML, bfunc){
  size <- length(SBML$genes)
  if (bfunc > size | bfunc < 1){return(paste("Error: bfunc must be in range [1,", size, "].", sep=""))}
  #generate all relevant 2^input combinations 
  inputs <- SBML$interactions[[bfunc]]$input
  #print(paste("Function", bfunc, "has", length(inputs), "inputs."))
  allinputcombis <- generateuniquestates(N=2^length(inputs), l=length(inputs), returnList = F)
  colnames(allinputcombis) <- inputs
  #feed into bfunc, get mapping to 1 or 0
  #stringFunc <- SBML$interactions[[bfunc]]$expression
  outputs <- rep(NA, dim(allinputcombis)[1])
  #take allinputcombis[1,] -> map to vector of zeros except with indices "inputs" replaced with these values?
  for (i in 1:dim(allinputcombis)[1]){
    initialvec <- rep(0, size)
    initialvec[c(inputs)] <- allinputcombis[i,]
    initialvec <- unlist(initialvec)
    newstate <- unlist(stateTransition(SBML, state = initialvec, type="synchronous"))
    #Does bfunc(inputs==newstate) = 1?
    outputs[i] <- newstate[bfunc]
  }
  #delete all input combination from matrix which mapped to 0
  #remainder is supp f_bfunc (input1, input2,...)
  support <- allinputcombis[which(outputs == 1),]
  return(support)
}
sumsupp <- function(support, j, xj){
  if (xj %in% c(1,0) == F){return("Error: xj must be binary.")}
  #Get sum_x € supp f_i P(X=x | Xj=xj)
  #P(X=x | Xj=xj) = P(X=x, xj=xj)/P(xj=xj)
  Pxj <- 0.5 #by assumption, denominator of Bayesian expression always 1/2 (no bias)
  inputs <- dim(support)[2]
  #For how many B^n states in the support of bfunc is the value at position xj=xj?
  xjmatched <- length(which(support[, which(colnames(support) == j)] == xj))/(2^inputs)
  return(xjmatched/Pxj)
}
MutInfo <- function(SBML, fi, j){
  #MI(fi(X),Xj) = h(sumsuppfi p_x) - 
  # 0.5*h(sumsuppfi P(X=x|Xj=1)) - 
  # 0.5*h(sumsuppfi P(X=x|Xj=0))
  support <- suppBoolFunc(SBML, fi)
  if (typeof(support) == "integer"){ #what happens if BoolFunc fi has only one input, either pos. or neg.?
    term1 <- binaryShannonEntropy(1/2)
    #sum over suppfi px is 1/2 as there is one state in the support (either (1) or (0)) out of 2 possible states -> h(1/2) = 1
    term2 <- term3 <- binaryShannonEntropy(1)
    #either term2 (for neg. reg.) or term (for pos.reg.) don't exist anyway, only one element in support, only one input j
    #this input is either 1 or 0 once, the other case Xj=0 or Xj=1 never occurs
    #for the remaining term, e.g. term2: P(X=(1)|Xj=1) = P(1)/P(Xj=1) = 0.5/0.5 = 1 -> h(1)=0
    #--> single input BoolFuncs fi with only one regulator j always yield MutInfo=1-0.5*0-0.5*0=1
  } else {
    term1 <- binaryShannonEntropy(dim(support)[1]/2^(dim(support)[2]))
    #term1 = sum x € suppfi p_x = total nr of elements in suppfi / 2^inputs of fi
    term2 <- binaryShannonEntropy(sumsupp(support, j=j, xj=1))
    term3 <- binaryShannonEntropy(sumsupp(support, j=j, xj=0))
  }
  return(term1-0.5*term2-0.5*term3)
}
DP_j <- function(SBML, j){
  #get sum of MutInfo(fi(X), Xj) for all i that are regulated by j, i.e. all nodes having j as an input
  adjmat <- conv2adjmat(SBML)
  inputs <- which(adjmat[j,] == 1) #which nodes have j as a common input? need fi from these, given Xj=1 or Xj=0...
  #print(inputs)
  #print(paste("Function j=", j, " has inputs: ", toString(inputs), ".", sep=""))
  res <- 0
  for (i in inputs){
    MI_fi_Xj <- MutInfo(SBML, fi = i, j = j)
    res <- res +  MI_fi_Xj
  }
  return(res)
}
DeterminativePower <- function(SBML){
  size <- length(SBML$genes)
  DPvec <- rep(NA, size)
  for (j in 1:size){
    DPvec[j] <- DP_j(SBML, j)
  }
  names(DPvec) <- SBML$genes
  return(DPvec)
}

nroutlinks <- function(SBML){
  size <- length(SBML$genes)
  outlinks <- rep(NA, size)
  adjmat <- conv2adjmat(SBML)
  for (g in 1:size){
    outlinks[g] <- sum(adjmat[g,])
  }
  names(outlinks) <- SBML$genes
  return(outlinks)
}

#Generate z-scores (connectivity) normed by the network itself (not over all networks)
singlemoduleZ <- function(SBML){
  adjmat <- conv2adjmat(SBML)
  size <- length(SBML$genes)
  k <- rep(NA, size)
  z <- rep(NA, size)
  for (i in 1:size){
    k[i] <- sum(adjmat[i,]) + sum(adjmat[,i])
  }
  Kavg <- mean(k)
  Ksigma <- sd(k)
  
  for (i in 1:size){z[i] <- (k[i]-Kavg)/Ksigma}
  names(z) <- SBML$genes
  return(z)
}

# Comparison of scores in VBnDP and connectivity
percentile <- function(x, perc, fullrange=T){
  #fullrange=F => works as before but with sign changed! --> ecdf(-x)(-perc)
  #fullrange=T => highest value in x or above always yields 1, lowest value in x or below always yields 0 ==> NEW SCALING!
  if (fullrange==F){
    return(ecdf(-x)(-perc))
  }
  else {
    #get number of different values in vector x, find position # pos where perc falls in here -> 1 - (#pos - 1)/(nrunique-1)
    nrunique <- length(unique(x))
    if (nrunique == 1){
      print("Warning: No variation is given vector, only 1 unique score.")
      return(1)
      }
    sortedx <- sort(x)
      #perc <- 11.5 #for testing function
    if (perc %in% sortedx){
      nrpos <- which(unique(sortedx) == perc)
    }
    else {
      perc <- max(x[x <= perc]) #whatever actual entry in the vector is closest to the given value of perc (round down)
      nrpos <- which(unique(sortedx) == perc)
    }
    percscore <- 1 - (nrpos - 1)/(nrunique - 1)
    #print(percscore)
    return(percscore)
  }
}
