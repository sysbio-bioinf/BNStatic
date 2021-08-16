#Load necessary libraries, setwd #####
  rm(list=ls(all=TRUE))
  if(!require("BoolNet"))
  install.packages("BoolNet")
  if(!require("igraph"))
    install.packages("igraph")
  if(!require("ggplot2"))
    install.packages("ggplot2")
  if(!require("Matrix"))
    install.packages("Matrix")
  if(!require("ggpubr"))
    install.packages("ggpubr")
  if(!require("gtools"))
    install.packages("gtools")

#pathtoscripts <- paste0(dirname(rstudioapi::getSourceEditorContext()$path), "/") #Path to plotscripts + data folder
pathtoscripts <- "./"
source(paste0(pathtoscripts, "functions.R"))
path <- paste0(pathtoscripts, "Networks/")
nets <- mixedsort(dir(path, pattern = ".txt"))
#

##### Optional: generate data for static and dynamic measures for all networks ######
# Get nets, save nets in network measures in order to identify enumeration of nets
# Calculate measures and save in Networks/Network measures n_Z, n_VB, n_DP, n_Hamming, n_Loss, n_Gain
print(paste(length(nets), "networks found, calculating static and dynamic measures"))
saveRDS(nets, file=paste0(path, "/Network measures/nets.RDS")) #Save reference for indexing of nets

VBresults <- calculateAllMeasures(nets, path, measure2calc="VB", saveResults=TRUE, savepath=paste0(path, "/Network measures/"))
DPresults <- calculateAllMeasures(nets, path, measure2calc="DP", saveResults=TRUE, savepath=paste0(path, "/Network measures/"))
Zresults <- calculateAllMeasures(nets, path, measure2calc="Z", saveResults=TRUE, savepath=paste0(path, "/Network measures/"))
Hammingresults <- calculateAllMeasures(nets, path, measure2calc="Hamming", saveResults=TRUE, savepath=paste0(path, "/Network measures/"))
AttrLossresults <- calculateAllMeasures(nets, path, measure2calc="AttrLoss", saveResults=TRUE, savepath=paste0(path, "/Network measures/"))
AttrGainresults <- calculateAllMeasures(nets, path, measure2calc="AttrGain", saveResults=TRUE, savepath=paste0(path, "/Network measures/"))

##### Optional: Load results of static and dynamic measures, determines selection threshold for VBnDP #####
#Vary threshold T from 1 to 100, get value of T where sensitivity and specificity intersect as well as their value at this point
measure <- "VBnDP" #intersection
plotAvgDynSensSpec(nets, path, measure = measure, loaddata=TRUE)

##### Return suggested intervention targets #####
#Use previously determined threshold T for selecting nodes in VBnDP
#For the chosen network, returns a list of PM nodes as possible intervention targets
#Nodes are labelled as Hubs or Non-Hubs and sorted from the largest to smallest mismatch between VBnDP and connectivity
  nets
  threshold <- 73
  # Network index n, return target suggestions for this network
  suggestTargets(n=21, threshold=threshold, loaddata=TRUE)
  
  