#Load necessary libraries, setwd #####
rm(list=ls(all=TRUE))
if(!require("BoolNet"))
  install.packages("BoolNet")
library("BoolNet")
library("igraph")
library("ggplot2")
library("Matrix")
library("ggpubr")
library("gtools")

pathtoscripts <- utils::getSrcDirectory #Path to plotscripts + data folder
source(paste0(pathtoscripts, "functions.R"))
path <- paste0(pathtoscripts, "Networks/")
nets <- mixedsort(dir(path, pattern = ".txt"))
#

##### Calculate static and dynamic measures for all networks, save data #####
# Generate results of static and dynamic measures across networks, alternative load from results folder
for (n in 1:length(nets)){
  print(paste0("Net nr. ", n))
  netfile <- paste(path, nets[n], sep = "")
  SBML <- loadNetwork(netfile)
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
  print("Calculating Hamming distance")
  Hamming <- meanofminattrdistsmerged(SBML, method="trinary") #Hamming distance
  saveRDS(Hamming, file=paste0(path, "/Network measures/", "net_", n, "_Hamming.RDS"))
  print("Calculating attractor loss")
  AttrLoss <- attrcoveragemerged(SBML, merging="min", method="trinary") #Attractor loss
  saveRDS(AttrLoss, file=paste0(path, "/Network measures/", "net_", n, "_AttrLoss.RDS"))
  print("Calculating attractor gain")
  AttrGain <- newattrgainbygenemerged(SBML, method="trinary") #Attractor gain
  saveRDS(AttrGain, file=paste0(path, "/Network measures/", "net_", n, "_AttrGain.RDS"))
}

##### Calculate sensitivity & specificity for every threshold of selection #####
#sensit, specif vectors contain sensitivity and specificity for every threshold T for the given network 
#in the comparison of the chosen static and dynamic measures
xs <- 4 #1=VB, 2=DP, 3=Union, 4=Intersection
ys <- 1 #1=Hamming distance, 2=Attractor loss, 3=Attractor gain
#
#Get sens/spec for t=1,...,100 for a given dynfunc when using only a single network -> for all networks
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
      Hamming <- readRDS(paste0(path, "/Network measures/", "net_", n, "_Hamming.RDS")) #Hamming distance
      AttrLoss <- readRDS(paste0(path, "/Network measures/", "net_", n, "_AttrLoss.RDS")) #Attractor loss
      AttrGain <- readRDS(paste0(path, "/Network measures/", "net_", n, "_AttrGain.RDS")) #Attractor gain
      
      #if picking given percentage of high scoring genes instead of mean or mean+sd as threshold
      p <- ceiling((toppercent/100)*dim(adjmat)[1])
      
      #VB
      x <- VB
      xmean <- mean(x)
      xsd <- sd(x)
      xselect_s4 <- names(sort(x, decreasing = TRUE)[1:p])
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
      if (xs == 1){xselect <- xselect_s4}  #VB
      if (xs == 2){xselect <- xselect_DP}  #DP
      if (xs == 3){xselect <- union(xselect_s4, xselect_DP)}  #VB v DP
      if (xs == 4){xselect <- intersect(xselect_s4, xselect_DP)}  #VB n DP
      
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

##### Calculate max. MI along all simple paths to Hubs (Fig3D) #####
# Pick a non-hub starting node and a hub end node. Calculate all simple paths connecting these. Average the Mutual Information
# along the edges in these paths, find the path with the maximal MI. Average the maximal MIs for all possible start/end node 
# combinations, normalise by the number of hubs and number of non-hubs of the chosen class in the network.
# ==> PM contain higher MI (i.e. "more determining") paths towards hubs than NM or NS non-hubs
HubDefinition <- "Guimera"

#Calculate all simple paths from non-hub nodes to hubs, get mutual information along these paths
#get maximal MI along paths starting and ending at the same nodes
toppercent <- threshold <- 73
maxMIbynet_PM <- maxMIbynet_NM <- maxMIbynet_NS <- rep(NA, length(nets))
for (n in 1:length(nets)){
  print(n)
  ##### get percentiles, VBnDP rankings etc #####
  netfile <- paste(path, nets[n], sep = "")
  SBML <- loadNetwork(netfile)
  adjmat <- conv2adjmat(SBML, inputcorrected = F)
  graph <- graph_from_adjacency_matrix(adjmat)
  p <- ceiling((threshold/100)*dim(adjmat)[1])
  VB <- readRDS(paste0(path, "/Network measures/", "net_", n, "_VB.RDS"))
  DP <- readRDS(paste0(path, "/Network measures/", "net_", n, "_DP.RDS"))
  Z <- readRDS(paste0(path, "/Network measures/", "net_", n, "_Z.RDS"))
  Hamming <- readRDS(paste0(path, "/Network measures/", "net_", n, "_Hamming.RDS")) #Hamming distance
  AttrLoss <- readRDS(paste0(path, "/Network measures/", "net_", n, "_AttrLoss.RDS")) #Attractor loss
  AttrGain <- readRDS(paste0(path, "/Network measures/", "net_", n, "_AttrGain.RDS")) #Attractor gain
  AvgDynRank <- colMeans(rbind(rank(-Hamming, ties.method ="min"), rank(AttrLoss, ties.method="min"), rank(-AttrGain, ties.method="min")), na.rm=TRUE)
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
  
  ##### define sets to compare for given property #####
  nsindices <- which(colnames(adjmat) %in% setdiff(SBML$genes, names(VBnDP)))
  hubindices <- setdiff(which(Z > 2.5), nsindices)
  dispropindices <- which(colnames(adjmat) %in% setdiff(names(surprise), colnames(adjmat)[c(hubindices, nsindices)]))
  discardindices <- setdiff(seq(1:length(SBML$genes)), c(hubindices, nsindices, dispropindices))
  
  ##### Max DP along simple paths normed by path length from Node->Hub #####
  print(paste(length(hubindices), "hubs in NET", n))
  if (length(hubindices)>0){
    #PM
    if (length(dispropindices)>0){
      sumMI <- 0
      for (g in dispropindices){
        #print(paste0("Gene index ", g))
        paths <- all_simple_paths(graph, from=g, to=hubindices, mode="out")
        #get last element for every [[p]] paths entry, only retain those where this is index=h
        
        for (h in hubindices){
          if (g == h){next}
          lastelems <- which(sapply(paths, tail, 1) == h)
          pathstoh <- rep(list(NA), length(lastelems))
          #print(paste(length(lastelems), "paths to hub", h, "out of", length(paths), "paths to all hubs"))
          print("PM")
          if (length(pathstoh)> 0){
            
            for (el in 1:length(lastelems)){
              pathstoh[[el]] <- paths[[lastelems[el]]]
            }
            
            maxMIs <- rep(NA, length(pathstoh))
            print(paste(length(pathstoh), "paths from index", g, "to", h))
            #print(maxMIs)
            for (pa in 1:length(pathstoh)){ #p=list of gene names in path
              p <- pathstoh[[pa]]
              MIpath <- 0
              nodenames <- names(p)
              nrnodes <- length(nodenames)
              #which(colnames(adjmat) %in% nodenames)
              nodenames
              nodeindices <- rep(NA, length(nodenames))
              for (v in 1:length(nodeindices)){
                nodeindices[v] <- which(colnames(adjmat) == nodenames[v])
              }# edge indices going from g->h
              for (e in 1:(length(nodeindices)-1)){
                MI <- MutInfo(SBML, fi=nodeindices[e+1], j=nodeindices[e])
                #sum up along path, avg by pathlength, then add to sumMI
                MIpath <- MIpath + MI
              }
              MIpath <- MIpath/(length(nodeindices)-1)
              maxMIs[pa] <- MIpath
              #print(MIpath)
            }#paths loop
            
            #print(maxMIs)
            
            sumMI <- sumMI + max(maxMIs)
          }#if length paths>0
        }#end hub loop
      }#end PM loop
      sumMI <- sumMI/(length(dispropindices)*length(hubindices))
      maxMIbynet_PM[n] <- sumMI
    }#if PMs exist
    
    #NM
    if (length(discardindices)>0){
      sumMI <- 0 
      for (g in discardindices){
        print(paste0("Gene index ", g))
        paths <- all_simple_paths(graph, from=g, to=hubindices, mode="out")
        #get last element for every [[p]] paths entry, only retain those where this is index=h
        
        for (h in hubindices){
          if (g == h){next}
          lastelems <- which(sapply(paths, tail, 1) == h)
          pathstoh <- rep(list(NA), length(lastelems))
          print(paste(length(lastelems), "paths to hub", h, "out of", length(paths), "paths to all hubs"))
          print("NM")
          if (length(pathstoh)> 0){
            
            for (el in 1:length(lastelems)){
              pathstoh[[el]] <- paths[[lastelems[el]]]
            }
            
            maxMIs <- rep(NA, length(pathstoh))
            #print(paste(length(pathstoh), "paths from index", g, "to", h))
            #print(maxMIs)
            for (pa in 1:length(pathstoh)){ #p=list of gene names in path
              p <- pathstoh[[pa]]
              MIpath <- 0
              nodenames <- names(p)
              nrnodes <- length(nodenames)
              #which(colnames(adjmat) %in% nodenames)
              nodenames
              nodeindices <- rep(NA, length(nodenames))
              for (v in 1:length(nodeindices)){
                nodeindices[v] <- which(colnames(adjmat) == nodenames[v])
              }# edge indices flowing from g->h
              for (e in 1:(length(nodeindices)-1)){
                MI <- MutInfo(SBML, fi=nodeindices[e+1], j=nodeindices[e])
                #sum these up along the path, avg by pathlength, then add to sumMI
                MIpath <- MIpath + MI
              }
              MIpath <- MIpath/(length(nodeindices)-1)
              maxMIs[pa] <- MIpath
              #print(MIpath)
            }#paths loop
            #print(maxMIs)
            sumMI <- sumMI + max(maxMIs)
          }#if length paths>0
        }#end hub loop
      }#end NM loop
      sumMI <- sumMI/(length(discardindices)*length(hubindices))
      maxMIbynet_NM[n] <- sumMI
    }#if NMs exist
    
    #NS
    if (length(nsindices)>0){
      sumMI <- 0
      for (g in nsindices){
        print(paste0("Gene index ", g))
        paths <- all_simple_paths(graph, from=g, to=hubindices, mode="out")
        #get last element for every [[p]] paths entry, only retain those where this is index=h
        
        for (h in hubindices){
          if (g == h){next}
          lastelems <- which(sapply(paths, tail, 1) == h)
          pathstoh <- rep(list(NA), length(lastelems))
          print(paste(length(lastelems), "paths to hub", h, "out of", length(paths), "paths to all hubs"))
          print("NS")
          if (length(pathstoh)> 0){
            
            for (el in 1:length(lastelems)){
              pathstoh[[el]] <- paths[[lastelems[el]]]
            }
            
            maxMIs <- rep(NA, length(pathstoh))
            #print(paste(length(pathstoh), "paths from index", g, "to", h))
            #print(maxMIs)
            for (pa in 1:length(pathstoh)){ #p=list of gene names in path
              p <- pathstoh[[pa]]
              MIpath <- 0
              nodenames <- names(p)
              nrnodes <- length(nodenames)
              #which(colnames(adjmat) %in% nodenames)
              nodenames
              nodeindices <- rep(NA, length(nodenames))
              for (v in 1:length(nodeindices)){
                nodeindices[v] <- which(colnames(adjmat) == nodenames[v])
              }# edge indices flowing from g->h
              for (e in 1:(length(nodeindices)-1)){
                MI <- MutInfo(SBML, fi=nodeindices[e+1], j=nodeindices[e])
                #sum these up along the path, avg by pathlength, then add to sumMI
                MIpath <- MIpath + MI
              }
              MIpath <- MIpath/(length(nodeindices)-1)
              maxMIs[pa] <- MIpath
              #print(MIpath)
            }#paths loop
            #print(maxMIs)
            sumMI <- sumMI + max(maxMIs)
          }#if length paths>0
        }#end hub loop
      }#end NS loop
      sumMI <- sumMI/(length(nsindices)*length(hubindices))
      maxMIbynet_NS[n] <- sumMI
    }#if NS exist
    
  }#end over all if hubs exist
  #
}
print(maxMIbynet_PM)
print(maxMIbynet_NM)
print(maxMIbynet_NS)


##### Plot Figure 2A/S3 #####
# Find selection threshold based on sensitivity/specificity of a static measure vs the average of three sens./spec. - curves
# measure <- "VB"
# measure <- "DP"
# measure <- "VBvDP" #union
measure <- "VBnDP" #intersection

if (measure == "VB"){
  titlestr <- "(A) Average of dynamic measures - VB"
  xs <- 1
} else if (measure == "DP"){
  titlestr <- "(B) Average of dynamic measures - DP"
  xs <- 2
} else if (measure == "VBvDP"){
  titlestr <- "(C) Average of dynamic measures - VB \u222A DP"
  xs <- 3
} else if (measure == "VBnDP"){
  titlestr <- "(D) Average of dynamic measures - VB \u2229 DP"
  xs <- 4
}

ysens_gain <- yspec_gain <- ysens_sd_gain <- yspec_sd_gain <- rep(0,100)
ysens_loss <- yspec_loss <- ysens_sd_loss <- yspec_sd_loss <- rep(0,100)
ysens_hd <- yspec_hd <- ysens_sd_hd <- yspec_sd_hd <- rep(0,100)
for (n in 1:length(nets)){
  #add up all gain*(1/length_nets)
  Hamming_sens <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/sensitivity_net", n, "_xs",xs, "_ys", 1 ,".RDS"))
  Hamming_spec <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/specificity_net", n, "_xs",xs, "_ys", 1 ,".RDS"))
  ysens_hd <- ysens_hd %+% (Hamming_sens*(1/length(nets)))
  yspec_hd <- yspec_hd %+% (Hamming_spec*(1/length(nets)))
  
  Loss_sens <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/sensitivity_net", n, "_xs",xs, "_ys", 2 ,".RDS"))
  Loss_spec <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/specificity_net", n, "_xs",xs, "_ys", 2 ,".RDS"))
  ysens_loss <- ysens_loss %+% (Loss_sens*(1/length(nets)))
  yspec_loss <- yspec_loss %+% (Loss_spec*(1/length(nets)))
  
  Gain_sens <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/sensitivity_net", n, "_xs",xs, "_ys", 3 ,".RDS"))
  Gain_spec <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/specificity_net", n, "_xs",xs, "_ys", 3 ,".RDS"))
  ysens_gain <- ysens_gain %+% (Gain_sens*(1/length(nets)))
  yspec_gain <- yspec_gain %+% (Gain_spec*(1/length(nets)))
}

ysens_sd_gain <- yspec_sd_gain <- ysens_sd_loss <- yspec_sd_loss <- ysens_sd_hd <- yspec_sd_hd <- rep(0,100)
ysens_gain_bynet <- yspec_gain_bynet <- ysens_loss_bynet <- yspec_loss_bynet <- ysens_hd_bynet <- yspec_hd_bynet <- rep(0, length(nets))
for (T in 1:100){
  print(T)
  for (n in 1:length(nets)){
    Hamming_sens <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/sensitivity_net", n, "_xs",xs, "_ys", 1 ,".RDS"))
    Hamming_spec <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/specificity_net", n, "_xs",xs, "_ys", 1 ,".RDS"))
    Loss_sens <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/sensitivity_net", n, "_xs",xs, "_ys", 2 ,".RDS"))
    Loss_spec <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/specificity_net", n, "_xs",xs, "_ys", 2 ,".RDS"))
    Gain_sens <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/sensitivity_net", n, "_xs",xs, "_ys", 3 ,".RDS"))
    Gain_spec <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/specificity_net", n, "_xs",xs, "_ys", 3 ,".RDS"))
  #loop over T: collect sensgains at every T, calc sd over all length_nets
    ysens_hd_bynet[n] <- Hamming_sens[T]
    yspec_hd_bynet[n] <- Hamming_spec[T]
    ysens_loss_bynet[n] <- Loss_sens[T]
    yspec_loss_bynet[n] <- Loss_spec[T]
    ysens_gain_bynet[n] <- Gain_sens[T]
    yspec_gain_bynet[n] <- Gain_spec[T]
    
    ysens_sd_hd[T] <- sd(ysens_hd_bynet, na.rm=T)
    yspec_sd_hd[T] <- sd(yspec_hd_bynet, na.rm=T)
    ysens_sd_loss[T] <- sd(ysens_loss_bynet, na.rm=T)
    yspec_sd_loss[T] <- sd(yspec_loss_bynet, na.rm=T)
    ysens_sd_gain[T] <- sd(ysens_gain_bynet, na.rm=T)
    yspec_sd_gain[T] <- sd(yspec_gain_bynet, na.rm=T)
  }
}


ysens <- rowMeans(cbind(ysens_gain, ysens_loss, ysens_hd))
yspec <- rowMeans(cbind(yspec_gain, yspec_loss, yspec_hd))
ysens_sd <- rowMeans(cbind(ysens_sd_gain, ysens_sd_loss, ysens_sd_hd))
yspec_sd <- rowMeans(cbind(yspec_sd_gain, yspec_sd_loss, yspec_sd_hd))

x <- seq(1:100)
df <- as.data.frame(cbind(x,ysens, yspec, ysens_sd, yspec_sd))
plot(
  ggplot(df, aes(x, y = ysens, color = "")) + 
    geom_point(aes(y = ysens, col = "Sensitivity", show.legend=F)) + labs(shape = "shape legend title", colour = "colour legend title") +
    geom_point(aes(y = yspec, col = "Specificity", show.legend=F)) + 
    geom_ribbon(data=df,aes(ymin=ysens-ysens_sd,ymax=ysens+ysens_sd, fill='SD (Sensitivity)'),alpha=0.2, color=NA, show.legend = F) + 
    geom_ribbon(data=df,aes(ymin=yspec-yspec_sd,ymax=yspec+yspec_sd, fill='SD (Specificity)'),alpha=0.2, color=NA, show.legend = F) + 
    scale_x_continuous(labels=c(0,"",20,"",40,"",60,"",80,"",100), breaks=c(0,10,20,30,40,50,60,70,80,90,100)) + 
    scale_y_continuous(labels=c(0,"",0.2,"",0.4,"",0.6,"",0.8,"",1), breaks=c(0,0.1,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1)) + 
    coord_cartesian(ylim = c(0, 1), xlim =c(0,100)) + theme(plot.title = element_text(hjust = 0.5)) +
    theme_grey(base_size = 25) + theme(legend.position = "none", legend.title=element_blank(),
                                       axis.text=element_text(size=25), axis.title=element_text(size=12.5), 
                                       plot.title = element_text(size=15)) +
    theme(panel.grid.major = element_line(colour="grey", size = (0.5)),
          panel.grid.minor = element_line(size = (0.0), colour="grey"), panel.background = element_blank()) + 
    labs(title=titlestr,
         x = expression(paste("Threshold percentage ", italic(T)," of genes labelled as 'high impact'")), 
         y = "Sensitivity & Specificity") +
    scale_colour_discrete("") + theme(plot.title = element_text(hjust = 0.5))
)

#

##### Plot Figure 2B/S4 #####
# Plot sensitivity-specificity curves using the comparison against individual dynamic measures
# measure <- "VB"
# measure <- "DP"
# measure <- "VBvDP" #union
measure <- "VBnDP" #intersection

# Load data: sensitivity and specificity for dynamic measures at every threshold
ysens_gain <- yspec_gain <- ysens_sd_gain <- yspec_sd_gain <- rep(0,100)
ysens_loss <- yspec_loss <- ysens_sd_loss <- yspec_sd_loss <- rep(0,100)
ysens_hd <- yspec_hd <- ysens_sd_hd <- yspec_sd_hd <- rep(0,100)
if (measure == "VB"){
  titlestr <- "(A) Individual dynamic measures - VB"
  xs <- 1
} else if (measure == "DP"){
  titlestr <- "(B) Individual dynamic measures - DP"
  xs <- 2
} else if (measure == "VBvDP"){
  titlestr <- "(C) Individual dynamic measures - VB \u222A DP"
  xs <- 3
} else if (measure == "VBnDP"){
  titlestr <- "(D) Individual dynamic measures - VB \u2229 DP"
  xs <- 4
}
for (n in 1:length(nets)){
  Hamming_sens <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/sensitivity_net", n, "_xs",xs, "_ys", 1 ,".RDS"))
  Hamming_spec <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/specificity_net", n, "_xs",xs, "_ys", 1 ,".RDS"))
  ysens_hd <- ysens_hd %+% (Hamming_sens*(1/length(nets)))
  yspec_hd <- yspec_hd %+% (Hamming_spec*(1/length(nets)))
  
  Loss_sens <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/sensitivity_net", n, "_xs",xs, "_ys", 2 ,".RDS"))
  Loss_spec <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/specificity_net", n, "_xs",xs, "_ys", 2 ,".RDS"))
  ysens_loss <- ysens_loss %+% (Loss_sens*(1/length(nets)))
  yspec_loss <- yspec_loss %+% (Loss_spec*(1/length(nets)))
  
  Gain_sens <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/sensitivity_net", n, "_xs",xs, "_ys", 3 ,".RDS"))
  Gain_spec <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/specificity_net", n, "_xs",xs, "_ys", 3 ,".RDS"))
  ysens_gain <- ysens_gain %+% (Gain_sens*(1/length(nets)))
  yspec_gain <- yspec_gain %+% (Gain_spec*(1/length(nets)))
}
ysens_gain
yspec_gain

ysens_sd_gain <- yspec_sd_gain <- ysens_sd_loss <- yspec_sd_loss <- ysens_sd_hd <- yspec_sd_hd <- rep(0,100)
ysens_gain_bynet <- yspec_gain_bynet <- ysens_loss_bynet <- yspec_loss_bynet <- ysens_hd_bynet <- yspec_hd_bynet <- rep(0, length(nets))
for (T in 1:100){
  print(T)
  for (n in 1:length(nets)){
    Hamming_sens <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/sensitivity_net", n, "_xs",xs, "_ys", 1 ,".RDS"))
    Hamming_spec <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/specificity_net", n, "_xs",xs, "_ys", 1 ,".RDS"))
    Loss_sens <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/sensitivity_net", n, "_xs",xs, "_ys", 2 ,".RDS"))
    Loss_spec <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/specificity_net", n, "_xs",xs, "_ys", 2 ,".RDS"))
    Gain_sens <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/sensitivity_net", n, "_xs",xs, "_ys", 3 ,".RDS"))
    Gain_spec <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/specificity_net", n, "_xs",xs, "_ys", 3 ,".RDS"))
    #loop over T: collect sensgains at every T, calc sd over all length_nets
    ysens_hd_bynet[n] <- Hamming_sens[T]
    yspec_hd_bynet[n] <- Hamming_spec[T]
    ysens_loss_bynet[n] <- Loss_sens[T]
    yspec_loss_bynet[n] <- Loss_spec[T]
    ysens_gain_bynet[n] <- Gain_sens[T]
    yspec_gain_bynet[n] <- Gain_spec[T]
    
    ysens_sd_hd[T] <- sd(ysens_hd_bynet, na.rm=T)
    yspec_sd_hd[T] <- sd(yspec_hd_bynet, na.rm=T)
    ysens_sd_loss[T] <- sd(ysens_loss_bynet, na.rm=T)
    yspec_sd_loss[T] <- sd(yspec_loss_bynet, na.rm=T)
    ysens_sd_gain[T] <- sd(ysens_gain_bynet, na.rm=T)
    yspec_sd_gain[T] <- sd(yspec_gain_bynet, na.rm=T)
  }
}

x <- seq(1:100)
df <- as.data.frame(cbind(x, ysens_hd, yspec_hd, ysens_loss, yspec_loss, ysens_gain, yspec_gain))

plot(
  ggplot(df, aes(x, y = ysens, color = "")) + 
    #gain
    geom_point(aes(y = ysens_gain, col = "Sensitivity", col="green"), col="#237d02", alpha=1, shape=19) +
    geom_point(aes(y = yspec_gain, col = "Specificity", col="green"), col="#237d02", alpha=1, shape=17) +
    #loss
    geom_point(aes(y = ysens_loss, col = "Sensitivity", col="orange"), col="#ab7105", alpha=1, shape=19) +
    geom_point(aes(y = yspec_loss, col = "Specificity", col="orange"), col="#ab7105", alpha=1, shape=17) +
    #Hamming
    geom_point(aes(y = ysens_hd, col = "Sensitivity", col="purple"), col="purple", alpha=1, shape=19) +
    geom_point(aes(y = yspec_hd, col = "Specificity", col="purple"), col="purple", alpha=1, shape=17) +
    
    scale_x_continuous(labels=c(0,"",20,"",40,"",60,"",80,"",100), breaks=c(0,10,20,30,40,50,60,70,80,90,100)) + 
    scale_y_continuous(labels=c(0,"",0.2,"",0.4,"",0.6,"",0.8,"",1), breaks=c(0,0.1,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1)) + 
    coord_cartesian(ylim = c(0, 1), xlim =c(0,100)) + theme(plot.title = element_text(hjust = 0.5)) +
    theme_grey(base_size = 25) + theme(legend.position = "bottom", 
                                       axis.text=element_text(size=25), axis.title=element_text(size=12.5), 
                                       plot.title = element_text(size=15)) + 
    theme(panel.grid.major = element_line(colour="grey", size = (0.5)),
          panel.grid.minor = element_line(size = (0.0), colour="grey"), panel.background = element_blank()) + 
    labs(title=titlestr,
         x = "Threshold percentage T of genes labelled as 'high impact'", y = "") + 
    scale_colour_discrete("") + theme(plot.title = element_text(hjust = 0.5))
)

#

##### Generate Plots for Figures 2A-D, comparison between classes #####
#All nodes are either a Hub or Non-Hub.
#Furthermore, all nodes are either not selected (NS), have positive mismatch (PM, VBnDP ranks higher than connectivity) 
# or no/negative mismatch (NM, VBnDP equal or lower ranking than connectivity)
#See lists of nodes according to both Hub definitions and their classifications in the given .xlsx tables 
# (Node classification tables)
#
##### Plot Figure 3A - Static impact #####
Figure <- "A"
HubDefinition <- "Guimera"
for (HubDefinition in c("Guimera")){
l <- length(nets)
toppercent <- 73
Havgperc <- nonHavgperc <- PMavgperc <- NMavgperc <- NSavgperc <- rep(NA, l)
Hcount <- nonHcount <- PMcount <- NMcount <- NScount <- 0

#Generate data
for (n in 1:length(nets)){
  print(n)
  ##### get measures #####
  # get percentiles, VBnDP rankings etc 
  netfile <- paste(path, nets[n], sep = "")
  SBML <- loadNetwork(netfile)
  adjmat <- conv2adjmat(SBML, inputcorrected = F)
  graph <- graph_from_adjacency_matrix(adjmat)
  size <- dim(adjmat)[1]
  p <- ceiling((toppercent/100)*size)
  VB <- readRDS(paste0(path, "Network measures/net_", n, "_VB.RDS"))
  DP <- readRDS(paste0(path, "Network measures/net_", n, "_DP.RDS"))
  Z <- readRDS(paste0(path, "Network measures/net_", n, "_Z.RDS"))
  Hamming <- readRDS(paste0(path, "Network measures/net_", n, "_Hamming.RDS")) #Hamming distance
  AttrLoss <- readRDS(paste0(path, "Network measures/net_", n, "_AttrLoss.RDS")) #Attractor loss
  AttrGain <- readRDS(paste0(path, "Network measures/net_", n, "_AttrGain.RDS")) #Attractor gain
  AvgDynRank <- colMeans(rbind(rank(-Hamming, ties.method ="min"), rank(AttrLoss, ties.method="min"), rank(-AttrGain, ties.method="min")), na.rm=TRUE)
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
    StatZQuantRatio[g] <- percentile(VBnDP, VBnDP[g], fullrange=T) - percentile(-ZmatchVBnDP, -ZmatchVBnDP[g], fullrange=T) 
    #low number=high rank
  }
  surprise <- which(StatZQuantRatio > 0)
  
  ##### define sets to compare for given property #####
  NSindices <- which(colnames(adjmat) %in% setdiff(names(AvgDynRank), names(VBnDP)))
  if (HubDefinition == "Guimera"){
    hubindices <- which(Z > 2.5) #Z>2.5 hub definition by Guimera et al.
  } else {print("Not a valid hub definition.")}
  nonhubindices <- setdiff(seq(1:length(SBML$genes)), hubindices)
  PMindices <- which(colnames(adjmat) %in% setdiff(names(surprise), colnames(adjmat)[c(NSindices)]))
  NMindices <- setdiff(seq(1:length(SBML$genes)), c(NSindices, PMindices))
  
  Hcount <- Hcount + length(hubindices)
  nonHcount <- nonHcount + length(nonhubindices)
  PMcount <- PMcount + length(PMindices)
  NMcount <- NMcount + length(NMindices)
  NScount <- NScount + length(NSindices)
  
  ##### write to file #####
  if (Figure == "A"){
    resvec <- sort(VBnDP)
  } else if (Figure %in% c("B","C")){
    resvec <- sort(AvgDynRank)
  }
  labelvec <- rep(NA, length(resvec))
  #get average percentile position of Hs, pHs, nphs in resvec -> print to three avgvecs
  Havgs <- nonHavgs <- PMavgs <- NMavgs <- NSavgs <- c()
  for (g in 1:length(resvec)){
    index <- which(colnames(adjmat) == names(resvec)[g])
    if (index %in% hubindices){
      labelvec[g] <- "Hub"
      if (Figure %in% c("A", "B")){
        Havgs <- append(Havgs, percentile(resvec, resvec[g]))
      } else if (Figure == "C"){
        Havgs <- append(Havgs, Z[index])
      }
    }#if in hubindices
    if (index %in% nonhubindices){
      labelvec[g] <- "Non-Hub"
      if (Figure %in% c("A", "B")){
        nonHavgs <- append(nonHavgs, percentile(resvec, resvec[g]))
      } else if (Figure == "C"){
        nonHavgs <- append(nonHavgs, Z[index])
      }
    }#if in nonhubindices
    if (index %in% PMindices){
      labelvec[g] <- "PM"
      if (Figure %in% c("A", "B")){
        PMavgs <- append(PMavgs, percentile(resvec, resvec[g]))
      } else if (Figure == "C"){
        PMavgs <- append(PMavgs, Z[index])
      }
    }
    if (index %in% NMindices){
      labelvec[g] <- "NM"
      if (Figure %in% c("A", "B")){
        NMavgs <- append(NMavgs, percentile(resvec, resvec[g]))
      } else if (Figure == "C"){
        NMavgs <- append(NMavgs, Z[index])
      }
    }
    if (index %in% NSindices){
      labelvec[g] <- "NS"
      if (Figure %in% c("A", "B")){
        NSavgs <- append(NSavgs, percentile(resvec, resvec[g]))
      } else if (Figure == "C"){
        NSavgs <- append(NSavgs, Z[index])
      }
    }
  }
  Havgperc[n] <- mean(Havgs)
  nonHavgperc[n] <- mean(nonHavgs)
  PMavgperc[n] <- mean(PMavgs)
  NMavgperc[n] <- mean(NMavgs)
  NSavgperc[n] <- mean(NSavgs)
}

df <- matrix(0, nrow = l*4, ncol=2)
df[1:l,1] <- Havgperc
df[(l+1):(2*l),1] <- nonHavgperc
df[(2*l+1):(3*l),1] <- PMavgperc
df[(3*l+1):(4*l),1] <- NMavgperc
df[1:l,2] <- "Hub"
df[(l+1):(2*l),2] <- "Non-Hub"
df[(2*l+1):(3*l),2] <- "PM"
df[(3*l+1):(4*l),2] <- "NM"
df <- as.data.frame(df)
df
names(df) <- c("percentiles", "group")
df$percentiles <- as.numeric(as.character(df$percentiles))
df$group <- factor(df$group, levels = c("Hub", "Non-Hub", "PM", "NM"))

padjvals <- compare_means(percentiles ~ group, df, method="wilcox.test", paired=T, p.adjust.method = "bonferroni") #get p.adj
for (p in 1:length(padjvals$p.adj)){
  if (padjvals$p.adj[p] == 1){padjvals$p.adj[p] <- ">0.99"}
}

if (HubDefinition == "Guimera"){
  pvals_Guimera <- padjvals
  Hcount_Guimera <- Hcount
  nonHcount_Guimera <- nonHcount
  PMcount_Guimera <- PMcount
  NMcount_Guimera <- NMcount
  NScount_Guimera <- NScount
  Havgperc_Guimera <- Havgperc
  nonHavgperc_Guimera <- nonHavgperc
  PMavgperc_Guimera <- PMavgperc
  NMavgperc_Guimera <- NMavgperc
  NSavgperc_Guimera <- NSavgperc
}

plot(
  ggboxplot(df, x = "group" , y = "percentiles", xlab="Classification", ylab="Percentile score in static impact ranking") + 
    stat_pvalue_manual(padjvals, label = "p.adj", y.position = c(1.1,1.2,1.3,
                                                                 1.5,1.6,
                                                                 1.8)) + ggtitle("   A") + 
    theme(plot.title = element_text(size = 20, face = "bold")) + 
    annotate("text", x=c(1,2,3,4), y=-0.05, 
             label=c(paste0("n=", Hcount),paste0("n=", nonHcount),paste0("n=", PMcount),paste0("n=", NMcount)))
)
}
#

##### Plot Figure 3B - Dynamic impact #####
# In the Guimera hub definition, all 17 hubs are in the NM class, giving this class higher dynamic impact,
# even though hubs have much higher connectivity than non-hub NM nodes
Figure <- "B"
HubDefinition <- "Guimera"
for (HubDefinition in c("Guimera")){
l <- length(nets)
toppercent <- 73
Havgperc <- nonHavgperc <- PMavgperc <- NMavgperc <- NSavgperc <- rep(NA, l)
Hcount <- nonHcount <- PMcount <- NMcount <- NScount <- 0


#Generate data
for (n in 1:length(nets)){
  print(n)
  ##### get measures #####
  # get percentiles, VBnDP rankings etc 
  netfile <- paste(path, nets[n], sep = "")
  SBML <- loadNetwork(netfile)
  adjmat <- conv2adjmat(SBML, inputcorrected = F)
  graph <- graph_from_adjacency_matrix(adjmat)
  size <- dim(adjmat)[1]
  p <- ceiling((toppercent/100)*size)
  VB <- readRDS(paste0(path, "Network measures/net_", n, "_VB.RDS"))
  DP <- readRDS(paste0(path, "Network measures/net_", n, "_DP.RDS"))
  Z <- readRDS(paste0(path, "Network measures/net_", n, "_Z.RDS"))
  Hamming <- readRDS(paste0(path, "Network measures/net_", n, "_Hamming.RDS")) #Hamming distance
  AttrLoss <- readRDS(paste0(path, "Network measures/net_", n, "_AttrLoss.RDS")) #Attractor loss
  AttrGain <- readRDS(paste0(path, "Network measures/net_", n, "_AttrGain.RDS")) #Attractor gain
  AvgDynRank <- colMeans(rbind(rank(-Hamming, ties.method ="min"), rank(AttrLoss, ties.method="min"), rank(-AttrGain, ties.method="min")), na.rm=TRUE)
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
    StatZQuantRatio[g] <- percentile(VBnDP, VBnDP[g], fullrange=T) - percentile(-ZmatchVBnDP, -ZmatchVBnDP[g], fullrange=T) 
    #low number=high rank
  }
  surprise <- which(StatZQuantRatio > 0)
  
  ##### define sets to compare for given property #####
  NSindices <- which(colnames(adjmat) %in% setdiff(names(AvgDynRank), names(VBnDP)))
  if (HubDefinition == "Guimera"){
    hubindices <- which(Z > 2.5) #Z>2.5 hub definition by Guimera et al.
  } else {print("Not a valid hub definition.")}
  nonhubindices <- setdiff(seq(1:length(SBML$genes)), hubindices)
  PMindices <- which(colnames(adjmat) %in% setdiff(names(surprise), colnames(adjmat)[c(NSindices)]))
  NMindices <- setdiff(seq(1:length(SBML$genes)), c(NSindices, PMindices))
  
  Hcount <- Hcount + length(hubindices)
  nonHcount <- nonHcount + length(nonhubindices)
  PMcount <- PMcount + length(PMindices)
  NMcount <- NMcount + length(NMindices)
  NScount <- NScount + length(NSindices)
  
  ##### write to file #####
  if (Figure == "A"){
    resvec <- sort(VBnDP)
  } else if (Figure %in% c("B","C")){
    resvec <- sort(AvgDynRank)
  }
  labelvec <- rep(NA, length(resvec))
  #get average percentile position of Hs, pHs, nphs in resvec -> print to three avgvecs
  Havgs <- nonHavgs <- PMavgs <- NMavgs <- NSavgs <- c()
  for (g in 1:length(resvec)){
    index <- which(colnames(adjmat) == names(resvec)[g])
    if (index %in% hubindices){
      labelvec[g] <- "Hub"
      if (Figure %in% c("A", "B")){
        Havgs <- append(Havgs, percentile(resvec, resvec[g]))
      } else if (Figure == "C"){
        Havgs <- append(Havgs, Z[index])
      }
    }#if in hubindices
    if (index %in% nonhubindices){
      labelvec[g] <- "Non-Hub"
      if (Figure %in% c("A", "B")){
        nonHavgs <- append(nonHavgs, percentile(resvec, resvec[g]))
      } else if (Figure == "C"){
        nonHavgs <- append(nonHavgs, Z[index])
      }
    }#if in nonhubindices
    if (index %in% PMindices){
      labelvec[g] <- "PM"
      if (Figure %in% c("A", "B")){
        PMavgs <- append(PMavgs, percentile(resvec, resvec[g]))
      } else if (Figure == "C"){
        PMavgs <- append(PMavgs, Z[index])
      }
    }
    if (index %in% NMindices){
      labelvec[g] <- "NM"
      if (Figure %in% c("A", "B")){
        NMavgs <- append(NMavgs, percentile(resvec, resvec[g]))
      } else if (Figure == "C"){
        NMavgs <- append(NMavgs, Z[index])
      }
    }
    if (index %in% NSindices){
      labelvec[g] <- "NS"
      if (Figure %in% c("A", "B")){
        NSavgs <- append(NSavgs, percentile(resvec, resvec[g]))
      } else if (Figure == "C"){
        NSavgs <- append(NSavgs, Z[index])
      }
    }
  }
  Havgperc[n] <- mean(Havgs)
  nonHavgperc[n] <- mean(nonHavgs)
  PMavgperc[n] <- mean(PMavgs)
  NMavgperc[n] <- mean(NMavgs)
  NSavgperc[n] <- mean(NSavgs)
}

#Fig3B
df <- matrix(0, nrow = l*5, ncol=2)
df[1:l,2] <- "Hub"
df[(l+1):(2*l),2] <- "Non-Hub"
df[(2*l+1):(3*l),2] <- "PM"
df[(3*l+1):(4*l),2] <- "NM"
df[(4*l+1):(5*l),2] <- "NS"
df[1:l,1] <- Havgperc
df[(l+1):(2*l),1] <- nonHavgperc
df[(2*l+1):(3*l),1] <- PMavgperc
df[(3*l+1):(4*l),1] <- NMavgperc
df[(4*l+1):(5*l),1] <- NSavgperc
df <- as.data.frame(df)
df
if (HubDefinition == "Guimera"){
  plabel <- c(1.15,1.3,1.5,1.7,
              1.9,2.1,2.3,
              2.5,2.7,2.9)
}
names(df) <- c("percentiles", "group")
df$percentiles <- as.numeric(as.character(df$percentiles))
df$group <- factor(df$group, levels = c("Hub", "Non-Hub", "PM", "NM", "NS"))
padjvals <- compare_means(percentiles ~ group, df, method="wilcox.test", paired=T, p.adjust.method = "bonferroni") #get p.adj
for (p in 1:length(padjvals$p.adj)){
  if (padjvals$p.adj[p] == 1){padjvals$p.adj[p] <- ">0.99"}
}

if (HubDefinition == "Guimera"){
  pvals_Guimera <- padjvals
  Hcount_Guimera <- Hcount
  nonHcount_Guimera <- nonHcount
  PMcount_Guimera <- PMcount
  NMcount_Guimera <- NMcount
  NScount_Guimera <- NScount
  Havgperc_Guimera <- Havgperc
  nonHavgperc_Guimera <- nonHavgperc
  PMavgperc_Guimera <- PMavgperc
  NMavgperc_Guimera <- NMavgperc
  NSavgperc_Guimera <- NSavgperc
}

plot(
  ggboxplot(df, x = "group" , y = "percentiles", xlab="Classification", ylab="Percentile score in dynamic impact ranking") + 
    stat_pvalue_manual(padjvals, label = "p.adj", y.position = plabel) + ggtitle("   B") + 
    theme(plot.title = element_text(size = 20, face = "bold")) + 
    annotate("text", x=c(1,2,3,4,5), y=-0, 
             label=c(paste0("n=", Hcount),paste0("n=", nonHcount),paste0("n=", PMcount),paste0("n=", NMcount),paste0("n=", NScount)))
)
}
#

##### Plot Figure 3C - Connectivity #####
Figure <- "C"
HubDefinition <- "Guimera"
for (HubDefinition in c("Guimera")){
l <- length(nets)
toppercent <- 73
Havgperc <- nonHavgperc <- PMavgperc <- NMavgperc <- NSavgperc <- rep(NA, l)
Hcount <- nonHcount <- PMcount <- NMcount <- NScount <- 0

#Generate data
for (n in 1:length(nets)){
  print(n)
  ##### get measures #####
  # get percentiles, VBnDP rankings etc 
  netfile <- paste(path, nets[n], sep = "")
  SBML <- loadNetwork(netfile)
  adjmat <- conv2adjmat(SBML, inputcorrected = F)
  graph <- graph_from_adjacency_matrix(adjmat)
  size <- dim(adjmat)[1]
  p <- ceiling((toppercent/100)*size)
  VB <- readRDS(paste0(path, "Network measures/net_", n, "_VB.RDS"))
  DP <- readRDS(paste0(path, "Network measures/net_", n, "_DP.RDS"))
  Z <- readRDS(paste0(path, "Network measures/net_", n, "_Z.RDS"))
  Hamming <- readRDS(paste0(path, "Network measures/net_", n, "_Hamming.RDS")) #Hamming distance
  AttrLoss <- readRDS(paste0(path, "Network measures/net_", n, "_AttrLoss.RDS")) #Attractor loss
  AttrGain <- readRDS(paste0(path, "Network measures/net_", n, "_AttrGain.RDS")) #Attractor gain
  AvgDynRank <- colMeans(rbind(rank(-Hamming, ties.method ="min"), rank(AttrLoss, ties.method="min"), rank(-AttrGain, ties.method="min")), na.rm=TRUE)
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
    StatZQuantRatio[g] <- percentile(VBnDP, VBnDP[g], fullrange=T) - percentile(-ZmatchVBnDP, -ZmatchVBnDP[g], fullrange=T) 
    #low number=high rank
  }
  surprise <- which(StatZQuantRatio > 0)
  
  ##### define sets to compare for given property #####
  NSindices <- which(colnames(adjmat) %in% setdiff(names(AvgDynRank), names(VBnDP)))
  if (HubDefinition == "Guimera"){
    hubindices <- which(Z > 2.5) #Z>2.5 hub definition by Guimera et al.
  } else {print("Not a valid hub definition.")}
  nonhubindices <- setdiff(seq(1:length(SBML$genes)), hubindices)
  PMindices <- which(colnames(adjmat) %in% setdiff(names(surprise), colnames(adjmat)[c(NSindices)]))
  NMindices <- setdiff(seq(1:length(SBML$genes)), c(NSindices, PMindices))
  
  Hcount <- Hcount + length(hubindices)
  nonHcount <- nonHcount + length(nonhubindices)
  PMcount <- PMcount + length(PMindices)
  NMcount <- NMcount + length(NMindices)
  NScount <- NScount + length(NSindices)
  
  ##### write to file #####
  if (Figure == "A"){
    resvec <- sort(VBnDP)
  } else if (Figure %in% c("B","C")){
    resvec <- sort(AvgDynRank)
  }
  labelvec <- rep(NA, length(resvec))
  #get average percentile position of Hs, pHs, nphs in resvec -> print to three avgvecs
  Havgs <- nonHavgs <- PMavgs <- NMavgs <- NSavgs <- c()
  for (g in 1:length(resvec)){
    index <- which(colnames(adjmat) == names(resvec)[g])
    if (index %in% hubindices){
      labelvec[g] <- "Hub"
      if (Figure %in% c("A", "B")){
        Havgs <- append(Havgs, percentile(resvec, resvec[g]))
      } else if (Figure == "C"){
        Havgs <- append(Havgs, Z[index])
      }
    }#if in hubindices
    if (index %in% nonhubindices){
      labelvec[g] <- "Non-Hub"
      if (Figure %in% c("A", "B")){
        nonHavgs <- append(nonHavgs, percentile(resvec, resvec[g]))
      } else if (Figure == "C"){
        nonHavgs <- append(nonHavgs, Z[index])
      }
    }#if in nonhubindices
    if (index %in% PMindices){
      labelvec[g] <- "PM"
      if (Figure %in% c("A", "B")){
        PMavgs <- append(PMavgs, percentile(resvec, resvec[g]))
      } else if (Figure == "C"){
        PMavgs <- append(PMavgs, Z[index])
      }
    }
    if (index %in% NMindices){
      labelvec[g] <- "NM"
      if (Figure %in% c("A", "B")){
        NMavgs <- append(NMavgs, percentile(resvec, resvec[g]))
      } else if (Figure == "C"){
        NMavgs <- append(NMavgs, Z[index])
      }
    }
    if (index %in% NSindices){
      labelvec[g] <- "NS"
      if (Figure %in% c("A", "B")){
        NSavgs <- append(NSavgs, percentile(resvec, resvec[g]))
      } else if (Figure == "C"){
        NSavgs <- append(NSavgs, Z[index])
      }
    }
  }
  Havgperc[n] <- mean(Havgs)
  nonHavgperc[n] <- mean(nonHavgs)
  PMavgperc[n] <- mean(PMavgs)
  NMavgperc[n] <- mean(NMavgs)
  NSavgperc[n] <- mean(NSavgs)
}
#

df <- matrix(0, nrow = l*5, ncol=2)
df[1:l,2] <- "Hub"
df[(l+1):(2*l),2] <- "Non-Hub"
df[(2*l+1):(3*l),2] <- "PM"
df[(3*l+1):(4*l),2] <- "NM"
df[(4*l+1):(5*l),2] <- "NS"
df[1:l,1] <- Havgperc
df[(l+1):(2*l),1] <- nonHavgperc
df[(2*l+1):(3*l),1] <- PMavgperc
df[(3*l+1):(4*l),1] <- NMavgperc
df[(4*l+1):(5*l),1] <- NSavgperc
df <- as.data.frame(df)
df
if (HubDefinition == "Guimera"){
  plabel <- c(4.2,4.6,5.0,5.4,
              1.4,1.8,2.3,
              2.7,3.1,3.7)
}

names(df) <- c("percentiles", "group")
df$percentiles <- as.numeric(as.character(df$percentiles))
df$group <- factor(df$group, levels = c("Hub", "Non-Hub", "PM", "NM", "NS"))
padjvals <- compare_means(percentiles ~ group, df, method="wilcox.test", paired=T, p.adjust.method = "bonferroni") #get p.adj
for (p in 1:length(padjvals$p.adj)){
  if (padjvals$p.adj[p] == 1){padjvals$p.adj[p] <- ">0.99"}
}

if (HubDefinition == "Guimera"){
  pvals_Guimera <- padjvals
  Hcount_Guimera <- Hcount
  nonHcount_Guimera <- nonHcount
  PMcount_Guimera <- PMcount
  NMcount_Guimera <- NMcount
  NScount_Guimera <- NScount
  Havgperc_Guimera <- Havgperc
  nonHavgperc_Guimera <- nonHavgperc
  PMavgperc_Guimera <- PMavgperc
  NMavgperc_Guimera <- NMavgperc
  NSavgperc_Guimera <- NSavgperc
}

plot(
  ggboxplot(df, x = "group" , y = "percentiles", xlab="Classification", ylab="Connectivity (z-score)") + 
    stat_pvalue_manual(padjvals, label = "p.adj", y.position = plabel) + ggtitle("   C") + 
    theme(plot.title = element_text(size = 20, face = "bold")) + 
    annotate("text", x=c(1,2,3,4,5), y=-2, 
             label=c(paste0("n=", Hcount),paste0("n=", nonHcount),paste0("n=", PMcount),paste0("n=", NMcount),paste0("n=", NScount)))
)
}
#


##### Plot Figure 3D - Max. MI along paths to Hubs #####
Figure <- "D"
HubDefinition <- "Guimera"
for (HubDefinition in c("Guimera", "Lu")){
l <- length(nets)
toppercent <- 73
Havgperc <- PMavgperc <- NMavgperc <- NSavgperc <- rep(NA, 35)
if (HubDefinition == "Guimera"){
  PMavgperc <- PMavgperc_Guimera <- readRDS(paste0(path, "Network measures/MaxMI/Guimera_maxMIbynet_PM.RDS"))
  NMavgperc <- NMavgperc_Guimera <- readRDS(paste0(path, "Network measures/MaxMI/Guimera_maxMIbynet_NM.RDS"))
  NSavgperc <- NSavgperc_Guimera <- readRDS(paste0(path, "Network measures/MaxMI/Guimera_maxMIbynet_NS.RDS"))
} else if (HubDefinition == "Lu"){
  PMavgperc <- PMavgperc_Lu <- readRDS(paste0(path, "Network measures/MaxMI/Lu_maxMIbynet_PM.RDS"))
  NMavgperc <- NMavgperc_Lu <- readRDS(paste0(path, "Network measures/MaxMI/Lu_maxMIbynet_NM.RDS"))
  NSavgperc <- NSavgperc_Lu <- readRDS(paste0(path, "Network measures/MaxMI/Lu_maxMIbynet_NS.RDS"))
}
df <- matrix(0, nrow = l*3, ncol=2)
df[1:l,2] <- "PM \n(Bottlenecks)"
df[(l+1):(2*l),2] <- "NM"
df[(2*l+1):(3*l),2] <- "NS \n(Non-selected)"
df[1:l,1] <- PMavgperc
df[(l+1):(2*l),1] <- NMavgperc
df[(2*l+1):(3*l),1] <- NSavgperc
df <- as.data.frame(df)
df
names(df) <- c("percentiles", "group")
df$percentiles <- as.numeric(as.character(df$percentiles))
df$group <- factor(df$group, levels = c("PM \n(Bottlenecks)", "NM", "NS \n(Non-selected)")) 

padjvals <- compare_means(percentiles ~ group, df, method="wilcox.test", paired=T, p.adjust.method = "bonferroni") #get p.adj

if (HubDefinition == "Guimera"){
  pvals_Guimera <- padjvals
} else if (HubDefinition == "Lu"){
  pvals_Lu <- padjvals
}
for (p in 1:length(padjvals$p.adj)){
  if (padjvals$p.adj[p] == 1){padjvals$p.adj[p] <- ">0.99"}
}

plot(
  ggboxplot(df, x = "group" , y = "percentiles", xlab="Classification", ylab="") + 
    stat_pvalue_manual(padjvals, label = "p.adj", y.position = c(0.7,0.8,0.9)) +
    ylab("Averaged maximal MI along paths to hubs") + 
    theme(plot.title = element_text(size = 20, face = "bold")) + ggtitle("   D") + 
    annotate("text", x=c(1,2,3), y=-0.1, 
             label=c(paste0("n=", PMcount),paste0("n=", NMcount),paste0("n=", NScount)))
)
}
#

##### Robustness tests for threshold #####
##### Bootstrap robustness test, plot FigS5 #####
#Loads seed for chosen random selection of 10000x35 networks, loads data for sensitivity & specificity, creates boxplot
#==> Selection threshold is robust against random selection of networks
# pick 35 nets out of 35 10000 times
# load sensitivity-specificity curves for given combinations of static and dynamic measures
# -> Determine thresholds each time
netseq <- seq(1:length(nets))
reps <- 10000

#If seed is already available for chosen networks, load here:
seed <- readRDS(paste0(pathtoscripts, "bootstrap/seed.RDS"))

#If seed does not yet exist, generate & save it
# seedmatrix <- matrix(NA, nrow = reps, ncol = length(nets))
# for (r in 1:reps){
#   randomselection <- sample(netseq, size=length(nets), replace = T)
#   seedmatrix[r,] <- randomselection
# }
# seed <- seedmatrix
# saveRDS(seedmatrix, file=paste0(pathtoscripts, "bootstrap/seed.RDS"))

nprime <- length(nets)
replacing = T #networks can be selected multiple times, nprime=n=35, replacing=T, select=F for bootstrap
select = F #== F if selecting, ==T if skipping the randomly chosen nets
gainT <- lossT <- hammingT <- balancedTs <- rep(NA,reps)
saving <- F #set TRUE if intermediate results should be saved
for (rand in 1:reps){
  print(paste0(rand, "/", reps))
  randvec <- seed[rand,]
  gain_avg_sens <- gain_avg_spec <- loss_avg_sens <- loss_avg_spec <- hamming_avg_sens <- hamming_avg_spec <- rep(0, 100)
  for (n in randvec){
    #ys: 1=Hamming, 2=AttrLoss, 3=AttrGain
    gain_netn_sens <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/sensitivity_net",n, "_xs4_ys3.RDS"))
    gain_netn_spec <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/specificity_net",n, "_xs4_ys3.RDS"))
    
    loss_netn_sens <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/sensitivity_net",n, "_xs4_ys2.RDS"))
    loss_netn_spec <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/specificity_net",n, "_xs4_ys2.RDS"))
    
    hamming_netn_sens <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/sensitivity_net",n, "_xs4_ys1.RDS"))
    hamming_netn_spec <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/specificity_net",n, "_xs4_ys1.RDS"))
    
    gain_avg_sens <- gain_avg_sens + gain_netn_sens
    gain_avg_spec <- gain_avg_spec %+% gain_netn_spec
    
    loss_avg_sens <- loss_avg_sens + loss_netn_sens
    loss_avg_spec <- loss_avg_spec %+% loss_netn_spec
    
    hamming_avg_sens <- hamming_avg_sens + hamming_netn_sens
    hamming_avg_spec <- hamming_avg_spec %+% hamming_netn_spec
    
  }
  gainT[rand] <- which(gain_avg_sens > gain_avg_spec)[1]
  lossT[rand] <- which(loss_avg_sens > loss_avg_spec)[1]
  hammingT[rand] <- which(hamming_avg_sens > hamming_avg_spec)[1]
  balancedTs[rand] <- mean(c(gainT[rand], lossT[rand], hammingT[rand]))
  if (saving == T | (rand==reps)){
  saveRDS(gainT, file=paste0(pathtoscripts, "bootstrap/gainT.RDS"))
  saveRDS(lossT, file=paste0(pathtoscripts, "bootstrap/lossT.RDS"))
  saveRDS(hammingT, file=paste0(pathtoscripts, "bootstrap/hammingT.RDS"))
  saveRDS(balancedTs, file=paste0(pathtoscripts, "bootstrap/balancedTs.RDS"))
  }
}

df <- as.matrix(c(gainT,lossT,hammingT,balancedTs))
groups <- c(rep("gainT", reps), rep("lossT", reps), rep("hammingT", reps), rep("balancedTs", reps))
df <- as.data.frame(cbind(df, groups))
names(df) <- c("Thresholds", "groups")

gainT <- readRDS(paste0(pathtoscripts, "bootstrap/gainT.RDS"))
lossT <- readRDS(paste0(pathtoscripts, "bootstrap/lossT.RDS"))
hammingT <- readRDS(paste0(pathtoscripts, "bootstrap/hammingT.RDS"))
balancedTs <- readRDS(paste0(pathtoscripts, "bootstrap/balancedTs.RDS"))
boxplot(gainT, lossT, hammingT, balancedTs, names=c("Attractor \ngain", "Attractor \nloss", "Hamming \ndistance", "Dynamic \naverage"),
        main="Bootstrap analysis of selection thresholds", cex.axis=0.8)
median(balancedTs)
IQR(balancedTs)
#

##### Generate data for permutation robustness test #####
#Check if there are multiple statically identical nodes at a given cutoff, test sensitivity/specificity for all combinations
#=> Selection threshold is robust in case of nodes with identical static ranking
xs <- 4 #1=VB, 2=DP, 3=union, 4=intersection
#ys <- 1 #1=Hamming distance, 2=Attractor loss, 3=Attractor gain
toppercent <- 73
#skipnets <- 14 #skip network if small variation of rankings leads to large number of combinations, save NA-vector instead
for (ys in 1:3){
  print(paste0("ys=",ys))
  sensit <- specif <- rep(NA, length(nets)) #write final averaged sens/spec values of a given net here
for (n in 1:length(nets)){
  sens100bynet <- spec100bynet <- rep(NA, 100)
  if (n %in% c(skipnets)){
    saveRDS(sens100bynet, file=paste0(path, "Network measures/SensitivitySpecificity/Figs4 data/sensitivity_net",n, "_xs",xs,"_ys",ys,".RDS"))
    saveRDS(spec100bynet, file=paste0(path, "Network measures/SensitivitySpecificity/Figs4 data/specificity_net",n, "_xs",xs,"_ys",ys,".RDS"))
    next
  }
  print(paste0("Net=", n))
  for (toppercent in 1:100){
    #print(paste0("Net=", n, ", p=",toppercent))
    netfile <- paste0(path, nets[n])
    SBML <- loadNetwork(netfile)
    adjmat <- conv2adjmat(SBML, inputcorrected = F)
    size <- dim(adjmat)[1]
    
    VB <- readRDS(paste0(path, "Network measures/net_", n, "_VB.RDS"))
    DP <- readRDS(paste0(path, "Network measures/net_", n, "_DP.RDS"))
    Hamming <- readRDS(paste0(path, "Network measures/net_", n, "_Hamming.RDS"))
    AttrLoss <- readRDS(paste0(path, "Network measures/net_", n, "_AttrLoss.RDS"))
    AttrGain <- readRDS(paste0(path, "Network measures/net_", n, "_AttrGain.RDS"))
    
    p <- ceiling((toppercent/100)*size) # p == T'
    
    xselect_VB <- sort(VB, decreasing = TRUE)[1:p]
    xselect_DP <- sort(DP, decreasing = TRUE)[1:p]
    yselect_Hamming <- sort(Hamming, decreasing = TRUE)[1:p]
    yselect_AttrLoss <- sort(AttrLoss)[1:p]
    yselect_AttrGain <- sort(AttrGain, decreasing = TRUE)[1:p]
    
    #What is the lowest selected value for all measures given this T'? ie. where is the cutoff?
    lowest_VB <- xselect_VB[length(xselect_VB)]
    lowest_DP <- xselect_DP[length(xselect_DP)]
    lowest_Hamming <- yselect_Hamming[length(yselect_Hamming)]
    lowest_AttrLoss <- yselect_AttrLoss[length(yselect_AttrLoss)]
    lowest_AttrGain <- yselect_AttrGain[length(yselect_AttrGain)]
    
    #how many entries are there in the full measures which have this exact value?
    #if it's more than 1, averages need to be taken for sens & spec at this specific T'
    nrcutVB <- length(which(VB == lowest_VB))
    nrcutDP <- length(which(DP == lowest_DP))
    nrcutHamming <- length(which(Hamming == lowest_Hamming))
    nrcutAttrLoss <- length(which(AttrLoss == lowest_AttrLoss))
    nrcutAttrGain <- length(which(AttrGain == lowest_AttrGain))
    
    #how many of these values are included given T' selections? i.e. choosing k out of n
    nrselVB <- length(which(xselect_VB == lowest_VB))
    nrselDP <- length(which(xselect_DP == lowest_DP))
    nrselHamming <- length(which(yselect_Hamming == lowest_Hamming))
    nrselAttrLoss <- length(which(yselect_AttrLoss == lowest_AttrLoss))
    nrselAttrGain <- length(which(yselect_AttrGain == lowest_AttrGain))
    
    #ys 1,2,3 == check overlap between VBnDP and d9,d13,d15 (Hamming, Loss, Gain) respectively
    if (ys == 1){
      dnum <- 9
      dynf <- Hamming
      lowest_y <- lowest_Hamming
      yselect_dynf <- yselect_Hamming
      nrseldyn <- nrselHamming
      totalDYNcombis <- choose(n = nrcutHamming, k = nrseldyn)
    } else if (ys == 2){
      dnum <- 13
      dynf <- AttrLoss
      lowest_y <- lowest_AttrLoss
      yselect_dynf <- yselect_AttrLoss
      nrseldyn <- nrselAttrLoss
      totalDYNcombis <- choose(n = nrcutAttrLoss, k = nrseldyn)
    } else if (ys == 3){
      dnum <- 15
      dynf <- AttrGain
      lowest_y <- lowest_AttrGain
      yselect_dynf <- yselect_AttrGain
      nrseldyn <- nrselAttrGain
      totalDYNcombis <- choose(n = nrcutAttrGain, k = nrseldyn)
    }
    dynmatrixcombis <- matrix(NA, nrow = totalDYNcombis, ncol=p)
    
    totalVBcombis <- choose(n = nrcutVB, k = nrselVB)
    totalDPcombis <- choose(n = nrcutDP, k = nrselDP)
    VBmatrixcombis <- matrix(NA, nrow = totalVBcombis, ncol=p)
    DPmatrixcombis <- matrix(NA, nrow = totalDPcombis, ncol=p)
    
    #Dynamic measures
    for (r in 1:dim(dynmatrixcombis)[1]){
      dynmatrixcombis[r,1:(p-nrseldyn)] <- names(yselect_dynf[1:(p-nrseldyn)])
    }
    borderlinegenes <- names(which(dynf == lowest_y)) 
    restcombs <- combinations(borderlinegenes, n = length(borderlinegenes), r = nrseldyn, repeats.allowed = F)
    dynmatrixcombis[,(p-nrseldyn+1):dim(dynmatrixcombis)[2]] <- restcombs
    
    sameVB <- xselect_VB[1:(p-nrselVB)]
    sameDP <- xselect_DP[1:(p-nrselDP)]
    
    #VB
    for (r in 1:dim(VBmatrixcombis)[1]){
      VBmatrixcombis[r,1:(p-nrselVB)] <- names(sameVB)
    }
    borderlinegenes <- names(which(VB == lowest_VB))
    restcombs <- combinations(borderlinegenes, n = length(borderlinegenes), r = nrselVB, repeats.allowed = F)
    VBmatrixcombis[,(p-nrselVB+1):dim(VBmatrixcombis)[2]] <- restcombs
    
    #DP
    for (r in 1:dim(DPmatrixcombis)[1]){
      DPmatrixcombis[r,1:(p-nrselDP)] <- names(sameDP)
    }
    borderlinegenes <- names(which(DP == lowest_DP)) 
    restcombs <- combinations(borderlinegenes, n = length(borderlinegenes), r = nrselDP, repeats.allowed = F)
    DPmatrixcombis[,(p-nrselDP+1):dim(DPmatrixcombis)[2]] <- restcombs
    
    maxintersectcombis <- dim(VBmatrixcombis)[1]*dim(DPmatrixcombis)[1]
    if (xs %in% c(1,2,4)){
      allintersectmatrix <- matrix(NA, nrow = 1, ncol=p)
    } else if (xs==3){
      allintersectmatrix <- matrix(NA, nrow = 1, ncol=length(SBML$genes))
    }
    
    for (i in 1:dim(VBmatrixcombis)[1]){
      for (j in 1:dim(DPmatrixcombis)[1]){
        if (xs==1){
          int <- VBmatrixcombis[i,]
        } else if (xs==2){
          int <- DPmatrixcombis[j,]
        } else if (xs==3){
          int <- union(VBmatrixcombis[i,], DPmatrixcombis[j,])
        } else if (xs==4){
          int <- intersect(VBmatrixcombis[i,], DPmatrixcombis[j,])
        }
        for (r in 1:dim(allintersectmatrix)[1]){
          identval <- setequal(int, allintersectmatrix[r,])
          if (identval == F & length(int > 0)){
            addvec <- rep(NA, dim(allintersectmatrix)[2])
            addvec[1:length(int)] <- int
            allintersectmatrix <- rbind(allintersectmatrix, addvec)
            break
          }
        }
      }
    }
    ### FORM ALL POSSIBLE OVERLAPS BETWEEN STATIC AND DYNAMIC MATRIX => Average sens & spec values ###
    if (dim(allintersectmatrix)[1]-1 > 0){
      #print(paste("Nr of VBnDP combis:", dim(allintersectmatrix)[1]-1))
      #print(paste("Nr of dynfunc combis:", dim(dynmatrixcombis)[1]))
      #print((dim(allintersectmatrix)[1]-1)*dim(dynmatrixcombis)[1])
      senslist <- speclist <- c() 
      if (dim(allintersectmatrix)[1] == 1){
        senslist <- append(senslist, NA)
        speclist <- append(speclist, NA)
      } else {
        for (r1 in 2:dim(allintersectmatrix)[1]){
          for (r2 in 1:dim(dynmatrixcombis)[1]){
            stat_noNA <- which(is.na(allintersectmatrix[r1,]) == F)
            s_high <- allintersectmatrix[r1,stat_noNA]
            s_low <- setdiff(colnames(adjmat), s_high)
            d_high <- dynmatrixcombis[r2,]
            d_low <- setdiff(colnames(adjmat), d_high)
            
            sensres <- length(intersect(s_high, d_high))/length(d_high)
            specres <- length(intersect(s_low, d_low))/length(d_low)
            
            if (is.na(sensres) == F){senslist <- append(senslist, sensres)}
            if (is.na(specres) == F){speclist <- append(speclist, specres)}
          }
          #  print(senslist)
          #  print(speclist)
        }
      }
    } else {senslist <- speclist <- c()}#end if VBnDP exists
    
    sens100bynet[toppercent] <- mean(senslist, na.rm=T)
    spec100bynet[toppercent] <- mean(speclist, na.rm=T)
    
    #save sens100bynet and spec100bynet vectors for every net, alternatively load from results folder
    
  }# end loop over T
  saveRDS(sens100bynet, file=paste0(path, "Network measures/SensitivitySpecificity/Figs4 data/sensitivity_net",n, "_xs",xs,"_ys",ys,".RDS"))
  saveRDS(spec100bynet, file=paste0(path, "Network measures/SensitivitySpecificity/Figs4 data/specificity_net",n, "_xs",xs,"_ys",ys,".RDS"))
}# #end loop over nets
}#end loop over dyn.function ys


### Plot FigS4
# measure <- "VB"
# measure <- "DP"
# measure <- "VBvDP" #union
measure <- "VBnDP" #intersection

# Load data: sensitivity and specificity for dynamic measures at every threshold
ysens_gain <- yspec_gain <- ysens_sd_gain <- yspec_sd_gain <- rep(0,100)
ysens_loss <- yspec_loss <- ysens_sd_loss <- yspec_sd_loss <- rep(0,100)
ysens_hd <- yspec_hd <- ysens_sd_hd <- yspec_sd_hd <- rep(0,100)
if (measure == "VB"){
  titlestr <- "(A) Individual dynamic measures - VB"
  xs <- 1
} else if (measure == "DP"){
  titlestr <- "(B) Individual dynamic measures - DP"
  xs <- 2
} else if (measure == "VBvDP"){
  titlestr <- "(C) Individual dynamic measures - VB \u222A DP"
  xs <- 3
} else if (measure == "VBnDP"){
  titlestr <- "(D) Individual dynamic measures - VB \u2229 DP"
  xs <- 4
}
for (n in 1:length(nets)){
  Hamming_sens <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/FigS4 data/sensitivity_net", n, "_xs",xs, "_ys", 1 ,".RDS"))
  Hamming_spec <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/FigS4 data/specificity_net", n, "_xs",xs, "_ys", 1 ,".RDS"))
  ysens_hd <- ysens_hd %+% (Hamming_sens*(1/length(nets)))
  yspec_hd <- yspec_hd %+% (Hamming_spec*(1/length(nets)))
  
  Loss_sens <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/FigS4 data/sensitivity_net", n, "_xs",xs, "_ys", 2 ,".RDS"))
  Loss_spec <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/FigS4 data/specificity_net", n, "_xs",xs, "_ys", 2 ,".RDS"))
  ysens_loss <- ysens_loss %+% (Loss_sens*(1/length(nets)))
  yspec_loss <- yspec_loss %+% (Loss_spec*(1/length(nets)))

  Gain_sens <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/FigS4 data/sensitivity_net", n, "_xs",xs, "_ys", 3 ,".RDS"))
  Gain_spec <- readRDS(paste0(path, "Network measures/SensitivitySpecificity/FigS4 data/specificity_net", n, "_xs",xs, "_ys", 3 ,".RDS"))
  ysens_gain <- ysens_gain %+% (Gain_sens*(1/length(nets)))
  yspec_gain <- yspec_gain %+% (Gain_spec*(1/length(nets)))
}

x <- seq(1:100)
df <- as.data.frame(cbind(x, ysens_hd, yspec_hd, ysens_loss, yspec_loss, ysens_gain, yspec_gain))
plot(
  ggplot(df, aes(x, y = ysens, color = "")) + 
    #gain
    geom_point(aes(y = ysens_gain, col = "Sensitivity", col="green"), col="#237d02", alpha=1, shape=19) +
    geom_point(aes(y = yspec_gain, col = "Specificity", col="green"), col="#237d02", alpha=1, shape=17) +
    #loss
    geom_point(aes(y = ysens_loss, col = "Sensitivity", col="orange"), col="#ab7105", alpha=1, shape=19) +
    geom_point(aes(y = yspec_loss, col = "Specificity", col="orange"), col="#ab7105", alpha=1, shape=17) +
    #Hamming
    geom_point(aes(y = ysens_hd, col = "Sensitivity", col="purple"), col="purple", alpha=1, shape=19) +
    geom_point(aes(y = yspec_hd, col = "Specificity", col="purple"), col="purple", alpha=1, shape=17) +
    
    scale_x_continuous(labels=c(0,"",20,"",40,"",60,"",80,"",100), breaks=c(0,10,20,30,40,50,60,70,80,90,100)) + 
    scale_y_continuous(labels=c(0,"",0.2,"",0.4,"",0.6,"",0.8,"",1), breaks=c(0,0.1,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1)) + 
    coord_cartesian(ylim = c(0, 1), xlim =c(0,100)) + theme(plot.title = element_text(hjust = 0.5)) +
    theme_grey(base_size = 25) + theme(legend.position = "bottom", 
                                       axis.text=element_text(size=25), axis.title=element_text(size=20), 
                                       plot.title = element_text(size=15)) + 
    theme(panel.grid.major = element_line(colour="grey", size = (0.5)),
          panel.grid.minor = element_line(size = (0.0), colour="grey"), panel.background = element_blank()) + 
    labs(title=titlestr,
         x = "Threshold percentage T of genes \nlabelled as 'high impact'", y = "") + 
    scale_colour_discrete("") + theme(plot.title = element_text(hjust = 0.5))
)

