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
library("xlsx")
library("ggVennDiagram")

pathtoscripts <- "./" #utils::getSrcDirectory() #Path to plotscripts + data folder
source(paste0(pathtoscripts, "functions.R"))
source(paste0(pathtoscripts, "pattern_templates.R"))
resultmatrix <- readRDS(paste0(pathtoscripts, "resultmatrix.RDS"))
path <- paste0(pathtoscripts, "Networks/")
nets <- mixedsort(dir(path, pattern = ".txt"))
#

##### Calculate static and dynamic measures for all networks, save data #####
# Generate results of static and dyna mic measures across networks, alternative load from results folder
# Save results of measures to savepath as RDS files if desired
VBresults <- calculateAllMeasures(nets, path, measure2calc="VB", saveResults=TRUE, savepath=paste0(path, "/Network measures/"))
DPresults <- calculateAllMeasures(nets, path, measure2calc="DP", saveResults=TRUE, savepath=paste0(path, "/Network measures/"))
Zresults <- calculateAllMeasures(nets, path, measure2calc="Z", saveResults=TRUE, savepath=paste0(path, "/Network measures/"))
Hammingresults <- calculateAllMeasures(nets, path, measure2calc="Hamming", saveResults=TRUE, savepath=paste0(path, "/Network measures/"))
AttrLossresults <- calculateAllMeasures(nets, path, measure2calc="AttrLoss", saveResults=TRUE, savepath=paste0(path, "/Network measures/"))
AttrGainresults <- calculateAllMeasures(nets, path, measure2calc="AttrGain", saveResults=TRUE, savepath=paste0(path, "/Network measures/"))
DynAvgresults <- calculateAllMeasures(nets, path, measure2calc="DynAvg", saveResults=TRUE, savepath=paste0(path, "/Network measures/"))

# Return top P percent of genes in a given network n according to some measure VB, DP, VBuDP, VBnDP, Connectivity, 
# Hamming, Loss, Gain, DynAvg
# Union and intersection of VB-DP can yield larger/smaller sets
n <- 3
P <- 73
returnTopRankings(measure="VBnDP", path, nets, n, P)

##### Calculate sensitivity & specificity for every threshold of selection #####
# sensit, specif vectors contain sensitivity and specificity for every threshold T for the given network 
# in the comparison of the chosen static and dynamic measures
# Get sens/spec for T=1,...,100 for a given dynfunc when using only a single network -> for all networks
HammingSensSpec <- calculateSensSpec(statmeasure="VBnDP", dynfunc="Hamming", nets, path, saveResults=TRUE, savepath=paste0(path, "/Network measures/"))
AttrLossSensSpec <- calculateSensSpec(statmeasure="VBnDP", dynfunc="AttrLoss", nets, path, saveResults=TRUE, savepath=paste0(path, "/Network measures/"))
AttrGainSensSpec <- calculateSensSpec(statmeasure="VBnDP", dynfunc="AttrGain", nets, path, saveResults=TRUE, savepath=paste0(path, "/Network measures/"))
# calculateSensSpec returns a list of length 2N for N networks in total
# The first N list entries contain sensitivities, the following N specificities
HammingSens <- HammingSensSpec[1:length(nets)]
HammingSpec <- HammingSensSpec[(length(nets)+1):(2*length(nets))]
AttrLossSens <- AttrLossSensSpec[1:length(nets)]
AttrLossSpec <- AttrLossSensSpec[(length(nets)+1):(2*length(nets))]
AttrGainSens <- AttrGainSensSpec[1:length(nets)]
AttrGainSpec <- AttrGainSensSpec[(length(nets)+1):(2*length(nets))]

##### Calculate max. MI along all simple paths to Hubs (Fig3D) #####
# Pick a non-hub starting node and a hub end node. Calculate all simple paths connecting these. Average the Mutual Information
# along the edges in these paths, find the path with the maximal MI. Average the maximal MIs for all possible start/end node 
# combinations, normalise by the number of hubs and number of non-hubs of the chosen class in the network.
# == > PM contain higher MI (i.e. "more determining") paths towards hubs than NM or NS non-hubs

# Calculate all simple paths from non-hub nodes to hubs, get mutual information along these paths
# Get maximal MI along paths starting and ending at the same nodes
toppercent <- threshold <- 73
maxMIs <- calcMaxMItoHubs(nets, path, threshold=threshold, loaddata=TRUE, saveaftereverynet = FALSE, saveResults = TRUE)
#


##### Plot Figure 2A #####
# Find selection threshold based on sensitivity/specificity of a static measure vs the average of three sens./spec. - curves
# measure <- "VB"
# measure <- "DP"
# measure <- "VBuDP" #union
measure <- "VBnDP" #intersection
plotAvgDynSensSpec(nets, path, measure = measure, loaddata=TRUE)
# 

##### Plot Figure 2B #####
# Plot sensitivity-specificity curves using the comparison against individual dynamic measures
# measure <- "VB"
# measure <- "DP"
# measure <- "VBuDP" #union
measure <- "VBnDP" #intersection
plotIndivDynSensSpec(nets, path, measure = measure, loaddata=TRUE)
#

##### Generate Plots for Figures 3A-D, comparison between classes #####
# All nodes are either a Hub or Non-Hub.
# Furthermore, all nodes are either not selected (NS), have positive mismatch (PM, VBnDP ranks higher than connectivity) 
# or no/negative mismatch (NM, VBnDP equal or lower ranking than connectivity)
# See lists of nodes according to both Hub definitions and their classifications in the given .xlsx table
# (Node classification tables)
#
##### Plot Figure 3A - Static impact #####
plotFig3ABC(nets, path, Figure="A", loaddata=TRUE)
#

##### Plot Figure 3B - Dynamic impact #####
# In the Guimera hub definition, all 21 hubs are in the NM class, giving this class higher dynamic impact,
# even though hubs have much higher connectivity than non-hub NM nodes
plotFig3ABC(nets, path, Figure="B", loaddata=TRUE)
#

##### Plot Figure 3C - Connectivity #####
plotFig3ABC(nets, path, Figure="C", loaddata=TRUE)
#

##### Plot Figure 3D - Max. MI along paths to Hubs #####
# Calculate paths first if loaddata=F
# If calculations have already been performed and saved as RDS, load these results instead with loaddata=T
plotFig3D(nets, path, loaddata=TRUE, threshold=73) 
#

##### Bootstrap robustness test, plot FigS5 #####
# Loads seed for chosen random selection of 10000x35 networks, loads data for sensitivity & specificity, creates boxplot
# ==> Selection threshold is robust against random selection of networks
# Pick 35 nets out of 35, 10000 times
# Load sensitivity-specificity curves for given combinations of static and dynamic measures
# -> Determine thresholds each time
plotThresholdBootstrap(nets, path, reps=10000, generateNewSeed=FALSE, loaddata=TRUE, saveResults=FALSE)
#
##### Generate feature matrix for machine learning algorithms #####
# Exemplary declaration of the first 100 nodes to be included in the test set, measures will be normalised according to the remaining training set
load("./Feature selection matrices/preparedDataSet_normedbynetsize.RData")
normaliseDataset(Dataset=dataset, test.ids = seq(1:100))
##### Plot Figure S6 - Overlap by network size #####
OverlapRatio <- sizes <- edges <- rep(NA, length(nets))
for (n in 1:length(nets)){
  print(paste0("NET: ", n))
  netfile <- paste(path, nets[n], sep = "")
  SBML <- loadNetwork(netfile)
  adjmat <- conv2adjmat(SBML, inputcorrected = F)
  sizes[n] <- dim(adjmat)[1]
  edges[n] <- sum(adjmat)
  TopStatOfNet <- returnTopRankings_readfiles(measure="VBnDP", path, nets, n, P = 73)
  TopDynOfNet <- returnTopRankings_readfiles(measure="DynAvg", path, nets, n, P = 73)
  matches <- intersect(names(TopStatOfNet), names(TopDynOfNet))
  unionStatDyn <- union(names(TopStatOfNet), names(TopDynOfNet))
  OverlapRatio[n] <- length(matches)/length(TopDynOfNet)
  
  print(OverlapRatio)
  print(mean(OverlapRatio, na.rm=TRUE))
}
#saveRDS(OverlapRatio, file=paste0(path, "StatDynOverlap_divby_StatSelectedSize.RDS"))

plot(x = sizes, y = OverlapRatio, xlab = "Size of network", ylab = "Overlap/Size of dynamic selection",
     main="Overlap between static and dynamic selection at T=73% as fraction of size of dynamic selection")
plot(x = (edges/sizes), y = OverlapRatio, xlab = "Number of edges/Size of network", ylab = "Overlap/Size of dynamic selection",
     main="Overlap between static and dynamic selection at T=73% as fraction of size of dynamic selection")


##### Plot Figure S7 - Network motif frequency by class #####
#resultmatrix.RDS contains classification results of all nodes across the 35 networks

GKboxCFFL <- NMboxCFFL <- NSboxCFFL <- rep(NA, length(nets))
GKboxIFFL <- NMboxIFFL <- NSboxIFFL <- rep(NA, length(nets))
GKboxBifan <- NMboxBifan <- NSboxBifan <- rep(NA, length(nets))

nodecounter <- 1
for (n in 1:length(nets)){
  print(n)
  netfile <- paste(path, nets[n], sep = "")
  netname <- sub('\\.txt$', '', nets[n]) 
  SBML <- loadNetwork(netfile)
  adjmat <- conv2adjmat(SBML)
  graph <- graph_from_adjacency_matrix(adjmat)
  trinadjmat <- trinaryadjmat(SBML)
  C1FFL_vec <- generalfflfinder(trinadjmat, ffltype = "C1")
  C2FFL_vec <- generalfflfinder(trinadjmat, ffltype = "C2")
  C3FFL_vec <- generalfflfinder(trinadjmat, ffltype = "C3")
  C4FFL_vec <- generalfflfinder(trinadjmat, ffltype = "C4")
  I1FFL_vec <- generalfflfinder(trinadjmat, ffltype = "I1")
  I2FFL_vec <- generalfflfinder(trinadjmat, ffltype = "I2")
  I3FFL_vec <- generalfflfinder(trinadjmat, ffltype = "I3")
  I4FFL_vec <- generalfflfinder(trinadjmat, ffltype = "I4")
  CFFL_vec <- C1FFL_vec + C2FFL_vec + C3FFL_vec + C4FFL_vec
  IFFL_vec <- I1FFL_vec + I2FFL_vec + I3FFL_vec + I4FFL_vec
  Bifan_vec <- motiffreqbygene(pattern_adjmat = bifan_adjmat, SBML = SBML)
  
  resultmatrix[6,nodecounter:(nodecounter+length(SBML$genes)-1)] <- CFFL_vec
  resultmatrix[7,nodecounter:(nodecounter+length(SBML$genes)-1)] <- IFFL_vec
  resultmatrix[8,nodecounter:(nodecounter+length(SBML$genes)-1)] <- Bifan_vec
  nodecounter <- nodecounter+length(SBML$genes)
  
  #Get corresponding class labels from resultmatrix[4,colinds]
  #which(CFFL_vec) are in GKindices, NMindices...
  colinds <- which(colnames(resultmatrix) == netname)
  classlabels <- resultmatrix[4,colinds]
  hublabels <- resultmatrix[5,colinds]
  indexGK <- which(classlabels == "Gatekeeper")
  indexNM <- which(classlabels == "NM")
  indexNS <- which(classlabels == "NS")
  indexHub <- which(hublabels == "Hub")
  indexNonHub <- which(hublabels == "Non-Hub")
  
  #Sum up participation in motif by class
  GKparticipationCFFL <- sum(CFFL_vec[indexGK]) 
  NMparticipationCFFL <- sum(CFFL_vec[indexNM])
  NSparticipationCFFL <- sum(CFFL_vec[indexNS])
  
  GKparticipationIFFL <- sum(IFFL_vec[indexGK])
  NMparticipationIFFL <- sum(IFFL_vec[indexNM])
  NSparticipationIFFL <- sum(IFFL_vec[indexNS])
  
  GKparticipationBifan <- sum(Bifan_vec[indexGK])
  NMparticipationBifan <- sum(Bifan_vec[indexNM])
  NSparticipationBifan <- sum(Bifan_vec[indexNS])
  
  #Results to plot, get motif participation for specific network
  GKboxCFFL[n] <- GKparticipationCFFL
  NMboxCFFL[n] <- NMparticipationCFFL
  NSboxCFFL[n] <- NSparticipationCFFL

  GKboxIFFL[n] <- GKparticipationIFFL
  NMboxIFFL[n] <- NMparticipationIFFL
  NSboxIFFL[n] <- NSparticipationIFFL

  GKboxBifan[n] <- GKparticipationBifan
  NMboxBifan[n] <- NMparticipationBifan
  NSboxBifan[n] <- NSparticipationBifan
}

mat <- matrix(NA, nrow = 35*9, ncol=4)
colnames(mat) <- c("Frequency", "Class", "Motif", "Selection")

mat[(0*35+1):(1*35),1] <- GKboxCFFL
mat[(1*35+1):(2*35),1] <- NMboxCFFL
mat[(2*35+1):(3*35),1] <- NSboxCFFL
mat[(0*35+1):(1*35),2] <- "PM \n (Gatekeeper)"
mat[(1*35+1):(2*35),2] <- "NM"
mat[(2*35+1):(3*35),2] <- "NS"
mat[1:(3*35),3] <- "C-FFL"
mat[(0*35+1):(1*35),4] <- "Selected"
mat[(1*35+1):(2*35),4] <- "Selected"
mat[(2*35+1):(3*35),4] <- "Non-selected"

mat[(3*35+1):(4*35),1] <- GKboxBifan
mat[(4*35+1):(5*35),1] <- NMboxBifan
mat[(5*35+1):(6*35),1] <- NSboxBifan
mat[(3*35+1):(4*35),2] <- "PM \n (Gatekeeper)"
mat[(4*35+1):(5*35),2] <- "NM"
mat[(5*35+1):(6*35),2] <- "NS"
mat[(3*35+1):(6*35),3] <- "Bifan"
mat[(3*35+1):(4*35),4] <- "Selected"
mat[(4*35+1):(5*35),4] <- "Selected"
mat[(5*35+1):(6*35),4] <- "Non-selected"

#Show occurrence of IFFL in same boxplot
mat[(6*35+1):(7*35),1] <- GKboxIFFL
mat[(7*35+1):(8*35),1] <- NMboxIFFL
mat[(8*35+1):(9*35),1] <- NSboxIFFL
mat[(6*35+1):(7*35),2] <- "PM \n (Gatekeeper)"
mat[(7*35+1):(8*35),2] <- "NM"
mat[(8*35+1):(9*35),2] <- "NS"
mat[(6*35+1):(9*35),3] <- "I-FFL"
mat[(6*35+1):(7*35),4] <- "Selected"
mat[(7*35+1):(8*35),4] <- "Selected"
mat[(8*35+1):(9*35),4] <- "Non-selected"

df <- as.data.frame(mat)
df[,1] <- as.numeric(df[,1])
df$Class <- factor(df$Class , levels=c("PM \n (Gatekeeper)", "NM", "NS"))
df$Selection <- factor(df$Selection , levels=c("Selected", "Non-selected"))

#PLOTTING:
#Selection (LOGSCALE)
ggplot(df,aes(x=Selection, y=Frequency, fill=Motif)) +
  geom_boxplot() + scale_y_log10()

#Class (LOGSCALE)
ggplot(df,aes(x=Class, y=Frequency, fill=Motif)) +
  geom_boxplot() + scale_y_log10()


##### Plot Figure S8 - Overlap of FVS with selection and classes #####
#Load netwise table containing classification and FVS: 
#Loop over nets -> check overlap of FVS by group, write to group specific vectors
netwiseTable <- readRDS(file=paste0(pathtoscripts, "/FVS/netwiseTable.RDS"))
Hub_FVS_overlapByNet <- PM_FVS_overlapByNet <- NM_FVS_overlapByNet <- NS_FVS_overlapByNet <- rep(NA,35)
FVSsizes <- rep(NA, 35)
for (n in 1:35){
  print(n)
  mat <- netwiseTable[[n]]
  Hubs <- na.omit(mat[,3])
  PMs <- na.omit(mat[,5])
  NMs <- na.omit(mat[,6])
  NSs <- na.omit(mat[,7])
  FVSs <- na.omit(mat[,8])
  FVSsizes[n] <- length(FVSs)
  if (length(Hubs) > 0){Hub_FVS_overlapByNet[n] <- length(intersect(FVSs, Hubs))/(length(FVSs))}
  if (length(PMs) > 0){PM_FVS_overlapByNet[n] <- length(intersect(FVSs, PMs))/(length(FVSs))}
  NM_FVS_overlapByNet[n] <- length(intersect(FVSs, NMs))/(length(FVSs))
  NS_FVS_overlapByNet[n] <- length(intersect(FVSs, NSs))/(length(FVSs))
}

#cols: Overlap, Class, Selection
mat <- matrix(NA, nrow=35*3, ncol=3)
colnames(mat) <- c("Overlap", "Class", "Selection")
mat[1:35,1] <- PM_FVS_overlapByNet
mat[(1*35+1):(2*35),1] <- NM_FVS_overlapByNet
mat[(2*35+1):(3*35),1] <- NS_FVS_overlapByNet
mat[1:35,2] <- "PM \n(Gatekeeper)"
mat[(1*35+1):(2*35),2] <- "NM"
mat[(2*35+1):(3*35),2] <- "NS"
mat[1:35,3] <- "Selected"
mat[(1*35+1):(2*35),3] <- "Selected"
mat[(2*35+1):(3*35),3] <- "Non-selected"
df <- as.data.frame(mat)
df[,1] <- as.numeric(df[,1])
df$Class <- factor(df$Class , levels=c("PM \n(Gatekeeper)", "NM", "NS")) #Correct order of boxes in plot
df$Selection <- factor(df$Selection , levels=c("Selected", "Non-selected")) #Correct order of boxes in plot

ggplot(df,aes(x=Selection, y=Overlap)) +
  geom_boxplot()

ggplot(df,aes(x=Class, y=Overlap)) +
  geom_boxplot()

#
##### Plot Figure S9 - Venn diagrams of Gatekeepers, motifs, canalyzers & FVS #####
#How often should a node participate in a motif/how many functions should it canalyze
CFFLthresh <- 1
IFFLthresh <- 1
Bifanthresh <- 1
Canalthresh <- 1

#Count FFLs and Bifan motif participation by gene across networks
totalIFFLs <- c()
for (n in 1:length(nets)){
  print(paste0("NET: ", n))
  netfile <- paste(path, nets[n], sep = "")
  netname <- sub('\\.txt$', '', nets[n]) 
  SBML <- loadNetwork(netfile)
  trinadjmat <- trinaryadjmat(SBML)
  #Sum up all four subtypes of coherent FFLs
  CFFL_vec <- generalfflfinder(trinadjmat, "C1") + generalfflfinder(trinadjmat, "C2") + 
    generalfflfinder(trinadjmat, "C3") + generalfflfinder(trinadjmat, "C4")
  #Sum up all four subtypes of incoherent FFLs
  IFFL_vec <- generalfflfinder(trinadjmat, "I1") + generalfflfinder(trinadjmat, "I2") + 
    generalfflfinder(trinadjmat, "I3") + generalfflfinder(trinadjmat, "I4")
  Bifan_vec <- motiffreqbygene(pattern_adjmat = bifan_adjmat, SBML = SBML)
}

#Alternatively, load results directly from resultmatrix
GKinds <- which(resultmatrix[4,] == "Gatekeeper")
CFFLinds <- which(as.numeric(resultmatrix[6,]) >= CFFLthresh)
IFFLinds <- which(as.numeric(resultmatrix[7,]) >= IFFLthresh)
Bifaninds <- which(as.numeric(resultmatrix[8,]) >= Bifanthresh)
Canalinds <- which(as.numeric(resultmatrix[9,]) >= Canalthresh)
FVSinds  <- which(as.numeric(resultmatrix[10,]) == 1) #1=FVS, 0=not in FVS
Motifinds <- union(union(CFFLinds, IFFLinds), Bifaninds)# nodes which participate in at least one motif of any kind

# Plot VennDiagram with overlap of GK, 3x motif, canalyzer, FVS (Figure S9, left)
x <- list(Gatekeeper=GKinds, Motif=Motifinds, Canalyzer=Canalinds, FVS=FVSinds)
venn <- Venn(x)
data <- process_data(venn)
ggVennDiagram(x, label = "count", label_size=5) + scale_fill_gradient(low="white",high = "red") +
  geom_sf(size=1, color = "grey", data = venn_setedge(data), show.legend = F)

# Plot VennDiagram with overlap of GK, individual motifs (Figure S9, right)
x <- list(Gatekeeper=GKinds, CFFL=CFFLinds, IFFL=IFFLinds, Bifan=Bifaninds)
venn <- Venn(x)
data <- process_data(venn)
ggVennDiagram(x, label = "count", label_size=5, set_size = 5) + scale_fill_gradient(low="white",high = "red") +
  geom_sf(size=1, color = "grey", data = venn_setedge(data), show.legend = F)


