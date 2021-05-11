#Load necessary libraries, setwd #####
rm(list=ls(all=TRUE))
library("BoolNet")
library("igraph")
library("ggplot2")
library("Matrix")
library("ggpubr")
library("gtools")

pathtoscripts <- paste0(dirname(rstudioapi::getSourceEditorContext()$path), "/") #Path to plotscripts + data folder
source(paste0(pathtoscripts, "functions.R"))
path <- paste0(pathtoscripts, "Networks/")
nets <- mixedsort(dir(path, pattern = ".txt"))
#

##### Optional: generate data for static and dynamic measures for all networks ######
#get nets, save nets in network measures in order to identify enumeration of nets
#calc measures and save in Networks/Network measures n_Z, n_VB, n_DP, n_Hamming, n_Loss, n_Gain
print(paste(length(nets), "networks found, calculating static and dynamic measures"))
saveRDS(nets, file=paste0(path, "/Network measures/nets.RDS")) #Save reference for indexing of nets

for (n in 1:length(nets)){
  print(paste0("Net nr. ", n))
  netfile <- paste(path, nets[n], sep = "")
  SBML <- loadNetwork(netfile) #need symbolic=TRUE if there are time delays
  adjmat <- conv2adjmat(SBML, inputcorrected = F)
  graph <- graph_from_adjacency_matrix(adjmat)
  #VB <- VertexBetweenness(params=list(adjmat=adjmat, normalized=T)) 
  print("Calculating static measures")
  VB <- VertexBetweenness(adjmat, normalized=T)
  saveRDS(VB, file=paste0(path, "/Network measures/", "net_", n, "_VB.RDS"))
  DP <- DeterminativePower(SBML)
  saveRDS(DP, file=paste0(path, "/Network measures/", "net_", n, "_DP.RDS"))
  Z <- singlemoduleZ(SBML)
  saveRDS(Z, file=paste0(path, "/Network measures/", "net_", n, "_Z.RDS"))
  #print("Calculating Hamming distance")
  #Hamming <- meanofminattrdistsmerged(SBML, method="trinary") #Hamming distance
  #saveRDS(Hamming, file=paste0(path, "/Network measures/", "net_", n, "_Hamming.RDS"))
  #print("Calculating attractor loss")
  #AttrLoss <- attrcoveragemerged(SBML, merging="min", method="trinary") #Attractor loss
  #saveRDS(AttrLoss, file=paste0(path, "/Network measures/", "net_", n, "_AttrLoss.RDS"))
  #print("Calculating attractor gain")
  #AttrGain <- newattrgainbygenemerged(SBML, method="trinary") #Attractor gain
  #saveRDS(AttrGain, file=paste0(path, "/Network measures/", "net_", n, "_AttrGain.RDS"))
}

##### Load results of static and dynamic measures, determine selection threshold for VBnDP #####

#Generate three sensitivity-specificity curves for VBnDP against Hamming distance, Attractor Loss, Attractor Gain across nets
# Save resulting sensitivity-specificity vectors for all nets
xs <- 4 #VBnDP as chosen static measure
length(nets)
for (n in 1:length(nets)){
  print(paste0("Net ", n))
  for (ys in 1:3){#Loop over all three dynamic measures
    print(paste0("ys=",ys))
      sensit <- rep(NA, 100)
      specif <- rep(NA, 100)
      for (toppercent in seq(1:100)){
        #print(paste("p=",toppercent,sep=""))
        netfile <- paste(path, nets[n], sep = "")
        SBML <- loadNetwork(netfile)
        adjmat <- conv2adjmat(SBML, inputcorrected = F)
        genenames <- colnames(adjmat)
        size <- dim(adjmat)[1]
        colnames(adjmat)
        
        VB <- readRDS(paste0(path, "/Network measures/", "net_", n, "_VB.RDS"))
        DP <- readRDS(paste0(path, "/Network measures/", "net_", n, "_DP.RDS"))
        Z <- readRDS(paste0(path, "/Network measures/", "net_", n, "_Z.RDS"))
        #Hamming <- readRDS(paste0(path, "/Network measures/", "net_", n, "_Hamming.RDS")) #Hamming distance
        #AttrLoss <- readRDS(paste0(path, "/Network measures/", "net_", n, "_AttrLoss.RDS")) #Attractor loss
        #AttrGain <- readRDS(paste0(path, "/Network measures/", "net_", n, "_AttrGain.RDS")) #Attractor gain
        
        #if picking given percentage of high scoring genes instead of mean or mean+sd as threshold
        p <- ceiling((toppercent/100)*dim(adjmat)[1])
        
        #VB
        x <- VB
        xmean <- mean(x)
        xsd <- sd(x)
        xselect_VB <- names(sort(x, decreasing = TRUE)[1:p])
        #DP
        x <- DP
        xmean <- mean(x)
        xsd <- sd(x)
        xselect_DP <- names(sort(x, decreasing = TRUE)[1:p])
        #Z
        x <- Z
        xmean <- mean(x)
        xsd <- sd(x)
        xselect_Z <- names(sort(x, decreasing = TRUE)[1:p])
        #Hamming
        y <- Hamming
        ymean <- mean(y)
        ysd <- sd(y)
        yselect_Hamming <- names(sort(y, decreasing = TRUE)[1:p])
        #AttrLoss
        y <- AttrLoss
        ymean <- mean(y)
        ysd <- sd(y)
        yselect_AttrLoss <- names(sort(y)[1:p])
        #AttrGain
        y <- AttrGain
        ymean <- mean(y)
        ysd <- sd(y)
        yselect_AttrGain <- names(sort(y, decreasing = TRUE)[1:p])
        
        HammingRank <- rank(-Hamming, ties.method = "min")
        AttrLossRank <- rank(AttrLoss, ties.method = "min")
        AttrGainRank <- rank(-AttrGain, ties.method = "min")
        
        ##### CHOOSE RELEVANT DYNAMIC MEASURE #####
        if (ys == 1){yselect <- yselect_Hamming}
        if (ys == 2){yselect <- yselect_AttrLoss}
        if (ys == 3){yselect <- yselect_AttrGain}
        ##### CHOOSE RELEVANT STATIC MESAURE #####
        if (xs == 1){xselect <- xselect_VB}  #VB
        if (xs == 2){xselect <- xselect_DP}  #DP
        if (xs == 3){xselect <- union(xselect_VB, xselect_DP)}  #VB v DP
        if (xs == 4){xselect <- intersect(xselect_VB, xselect_DP)}  #VB n DP
        
        s_high <- xselect
        s_low <- setdiff(genenames, xselect)
        d_high <- yselect
        d_low <- setdiff(genenames, yselect)
        
        sensit[toppercent] <- length(intersect(s_high, d_high))/length(d_high)
        specif[toppercent] <- length(intersect(s_low, d_low))/length(d_low)
        sensit[which(is.nan(sensit))] <- NA
        specif[which(is.nan(specif))] <- NA
        
      }#loop over threshold T
      #Option: Save generated vector sensit, specif
      saveRDS(sensit, file=paste0(path, "Network measures/SensitivitySpecificity/sensitivity_net", n, "_xs",xs, "_ys", ys ,".RDS"))
      saveRDS(specif, file=paste0(path, "Network measures/SensitivitySpecificity/specificity_net", n, "_xs",xs, "_ys", ys ,".RDS"))
  }#loop over dynamic functions
}#loop over nets
#

#Load curves for individual dynamic functions, average them
#return intersection point to be used as threshold for defining non-selected genes (NS)
xs <- 4 #VBnDP as static measure
HammingAvgSens <- HammingAvgSpec <- AttrLossAvgSens <- AttrLossAvgSpec <- AttrGainAvgSens <- AttrGainAvgSpec <- rep(0,100)
for (n in 1:length(nets)){
  #Load sens/spec for all nets for given dynfunc and average them
  HammingSens <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/sensitivity_net", n, "_xs",xs, "_ys", 1 ,".RDS"))
  HammingSpec <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/specificity_net", n, "_xs",xs, "_ys", 1 ,".RDS"))
  AttrLossSens <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/sensitivity_net", n, "_xs",xs, "_ys", 2 ,".RDS"))
  AttrLossSpec <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/specificity_net", n, "_xs",xs, "_ys", 2 ,".RDS"))
  AttrGainSens <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/sensitivity_net", n, "_xs",xs, "_ys", 3 ,".RDS"))
  AttrGainSpec <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/specificity_net", n, "_xs",xs, "_ys", 3 ,".RDS"))
  
  HammingAvgSens <- HammingAvgSens + HammingSens*(1/length(nets))
  HammingAvgSpec <- HammingAvgSpec %+% (HammingSpec*(1/length(nets)))
  AttrLossAvgSens <- AttrLossAvgSens + AttrLossSens*(1/length(nets))
  AttrLossAvgSpec <- AttrLossAvgSpec %+% (AttrLossSpec*(1/length(nets)))
  AttrGainAvgSens <- AttrGainAvgSens + AttrGainSens*(1/length(nets))
  AttrGainAvgSpec <- AttrGainAvgSpec %+% (AttrGainSpec*(1/length(nets)))
}
#Average the 3 curves of each dynamic measure, find sens>spec point, this yields selection threshold
AvgDynSens <- (HammingAvgSens+AttrLossAvgSens+AttrGainAvgSens)/3
AvgDynSpec <- (HammingAvgSpec+AttrLossAvgSpec+AttrGainAvgSpec)/3

SensGreaterThanSpec <- which(AvgDynSens > AvgDynSpec)
threshold <- SensGreaterThanSpec[1]
sensitivityatthreshold <- AvgDynSens[threshold]
print(paste0("The selection threshold for the static measure of VBnDP is determined to be at T=", threshold,
            "% with a sensitivity value of ",round(sensitivityatthreshold,3), " for the given set of ", length(nets), " networks."))

##### Choose Hub definition, return suggested intervention targets #####
HubDefinition <- "Guimera" #Hub = z-score > 2.5

threshold <- 73 #Use previously determined threshold T for selecting nodes in VBnDP
#For the chosen network, returns a list of PM nodes as possible intervention targets
#Nodes are labelled as Hubs or Non-Hubs and sorted from the largest to smallest mismatch between VBnDP and connectivity

#construct VBnDP and Z rankings for nodes in net, classify nodes
  nets
  n <-1 # Network index, return target suggestions for network nets[n]
  #get percentiles, rankings
  netfile <- paste(path, nets[n], sep = "")
  SBML <- loadNetwork(netfile)
  adjmat <- conv2adjmat(SBML, inputcorrected = F)
  p <- ceiling((threshold/100)*dim(adjmat)[1])
  VB <- readRDS(paste0(path, "/Network measures/", "net_", n, "_VB.RDS"))
  DP <- readRDS(paste0(path, "/Network measures/", "net_", n, "_DP.RDS"))
  Z <- readRDS(paste0(path, "/Network measures/", "net_", n, "_Z.RDS"))
  #Hamming <- readRDS(paste0(path, "/Network measures/", "net_", n, "_Hamming.RDS")) #Hamming distance
  #AttrLoss <- readRDS(paste0(path, "/Network measures/", "net_", n, "_AttrLoss.RDS")) #Attractor loss
  #AttrGain <- readRDS(paste0(path, "/Network measures/", "net_", n, "_AttrGain.RDS")) #Attractor gain
  #AvgDynRank <- colMeans(rbind(rank(-Hamming, ties.method ="min"), rank(AttrLoss, ties.method="min"), rank(-AttrGain, ties.method="min")), na.rm=TRUE)
  xselect_VB <- sort(VB, decreasing = TRUE)[1:p]
  xselect_DP <- sort(DP, decreasing = TRUE)[1:p]
  VBRank <- rank(-VB, ties.method="min")
  DPRank <- rank(-DP, ties.method="min")
  name <- intersect(names(xselect_VB), names(xselect_DP))
  VBnDP <- rep(NA, length(name))
  names(VBnDP) <- name
  Z_inVBnDP <- Z[which(names(Z) %in% names(VBnDP))]
  for (g in 1:length(VBnDP)){
    RankA <- VBRank[which(names(VBRank) == names(VBnDP)[g])]
    RankB <- DPRank[which(names(DPRank) == names(VBnDP)[g])]
    VBnDP[g] <- mean(c(RankA, RankB))
  }
  ZmatchVBnDP <- Z_inVBnDP[match(names(VBnDP), names(Z_inVBnDP))]
  StatZQuantRatio <- rep(NA, length(VBnDP))
  names(StatZQuantRatio) <- names(VBnDP)
  for (g in 1:length(VBnDP)){
    StatZQuantRatio[g] <- percentile(VBnDP, VBnDP[g]) - percentile(-ZmatchVBnDP, -ZmatchVBnDP[g]) #low number=high rank
  }
  surprise <- which(StatZQuantRatio > 0)
  
  #define sets
  NSindices <- which(colnames(adjmat) %in% setdiff(SBML$genes, names(VBnDP)))
  hubindices <- setdiff(which(Z > 2.5), NSindices)
  nonhubindices <- setdiff(seq(1:length(SBML$genes)), hubindices)
  PMindices <- which(colnames(adjmat) %in% setdiff(names(surprise), colnames(adjmat)[c(NSindices)]))
  NMindices <- setdiff(seq(1:length(SBML$genes)), c(NSindices, PMindices))

  targets <- colnames(adjmat)
  names(targets)[hubindices] <- "Hub"
  names(targets)[nonhubindices] <- "Non-Hub"
  targets <- targets[PMindices]
  decreasingMismatchOrder <- names(sort(StatZQuantRatio[surprise], decreasing = T))
  
  
  
  
  #Result: Recommended intervention targets, highest mismatch between VBnDP and connectivity rankings first
  targets <- targets[match(decreasingMismatchOrder, targets)]
  print(paste0("The recommended intervention targets (ordered by highest mismatch first) are:"))
  #Print with magnitude of mismatch:
  sort(StatZQuantRatio[surprise], decreasing = T)
  #Print with Hub/Non-Hub classification:
  print(targets)
  
  

  
  