setwd("Documents/Scolaire/ENS/Stage_M1/revbayes/tests/test_OBDP/output_inferKT_macro_test/")

library("ggplot2")
library("reshape2")
library("FossilSim")
library("RevGadgets")

S <- 200
weight_trees_posterior <- T

## Import Kt : probability distribution of the number of hidden lineages through time
Kt_trace_lines <- readLines("Kt_trace.txt")
Kt_trace_lines <- gsub("\\[| |\\]| ,|\t|;", "", Kt_trace_lines)[-1]               # Remove unwanted characters
Kt_trace <- read.csv(text = Kt_trace_lines, header = FALSE, na.strings = "nan")
Kt_trace[is.na(Kt_trace)] <- 0                                                    # Set NA values to 0
N <- length(Kt_trace)-1                                                           # Maximal number of hidden lineages
names(Kt_trace) <- 0:N                                                            # Set names to the number of hidden lineages

## Import the corresponding tree : get the number of observed lineages through time (LTT)
trees <- read.table("../output_macro_test_Mt/mcmc_OBDP_macro_test.trees", header = T)
trees$obd_tree <- sapply(trees$obd_tree, function(tree){read.tree(text=as.character(tree))})
burnin <- round(0.25*length(trees$Iteration))  # Change `0.25` value if done previously in RevBayes

## Add the iterations to Kt_trace
iterations <- trees$Iteration[(burnin+1):length(trees$Iteration)]
Kt_trace$Iteration <- rep(iterations, each=S)                                    # Add an iteration number column

## Browse all iterations
Kt_mean <- matrix(0, nrow=S, ncol=N+1)
observedLin_mean <- rep(0, S)
timePoints <- seq(0, -max(sapply(trees$obd_tree, tree.max)), length.out = S)
posteriors <- exp(trees$Posterior-max(trees$Posterior))
posteriors <- posteriors/sum(posteriors)
for (i in (burnin+1):length(trees$Iteration)){
  ### Increment the distribution of number of hidden lineades
  it <- trees$Iteration[i]
  print(it)
  Kt <- Kt_trace[Kt_trace$Iteration==it,-which(names(Kt_trace)=="Iteration")]
  lines_sum <- apply(Kt, 1, sum)
  for (j in 1:S){
    if (lines_sum[j]==0){
      Kt[j,1] <- 1.0                          # Empty lines are considered to have 0 hidden lineages
    }
    else if (lines_sum[j]!=0){
      Kt[j,] <- Kt[j,]/lines_sum[j]    # Normalise lines to 1 (several lines at 0.9999 or 1.0001)
    }
  }
  if (weight_trees_posterior){
    Kt_mean <- Kt_mean + Kt*posteriors[i]
  }
  else{ Kt_mean <- Kt_mean + Kt/length(iterations) }
  
  ### Get the LTT coordinates
  obd_tree <- collapse.singles(trees$obd_tree[[i]])    # Remove single nodes (ie. sampled ancestors)
  LTT <- data.frame(ltt.plot.coords(obd_tree))         # Extract LTT coordinates
  LTT$time <- round(LTT$time, 4)                       # Reduce precision (extant tips wrongly at time -0.000001)
  print(length(LTT$time))
  
  ### Increment number of observed lineages
  getNbObservedLin <- function(t){
    if (t < LTT$time[1]) return (0)
    first_consecutive_time_point <- which(LTT$time >= t)[1]
    return (LTT$N[first_consecutive_time_point])
  }
  observedLin <- sapply(timePoints, getNbObservedLin)
  if (weight_trees_posterior){
    observedLin_mean <- observedLin_mean + observedLin*posteriors[i]
  }
  else{ observedLin_mean <- observedLin_mean + observedLin/length(iterations) }
}

## Get the most probable number of hidden lineages (weighted mean according to their respective probabilities)
hiddenLin <- data.frame(weightedMean=matrix(apply(Kt_mean, 1, function(x){weighted.mean(as.integer(names(Kt_mean)), x)})),
                        maxProb=matrix(apply(Kt_mean, 1, function(x){as.integer(names(Kt_mean)[which(x==max(x))])})))
Kt_mean$aggregNbHiddenLin <- hiddenLin$weightedMean
# Kt_mean$aggregNbHiddenLin[is.na(Kt_mean$aggregNbHiddenLin)] <- 0  # Put NA values at 0 (`weighted.mean` artifact at t0)

## Get the aggregated number of observed and total lineages
Kt_mean$aggregNbTotalLin <- observedLin_mean + Kt_mean$aggregNbHiddenLin
Kt_mean$NbObservedLin <- round(observedLin_mean)

## Get the 95% credence interval
Kt_mean[1:51] <- t(apply(Kt_mean[1:51], 1, function(row){row/sum(row)}))   # Force lines to sum to 1 (correct numerical uncertainties)
Kt_mean_cumsum <- apply(Kt_mean[1:51], 1, cumsum)
Kt_mean$Cred0.025 <- observedLin_mean + apply(Kt_mean_cumsum, 2, function(col){which(col>0.025)[1]-1})
Kt_mean$Cred0.5 <- observedLin_mean + apply(Kt_mean_cumsum, 2, function(col){which(col>0.5)[1]})
Kt_mean$Cred0.975 <- observedLin_mean + apply(Kt_mean_cumsum, 2, function(col){which(col>0.975)[1]})

## Get the aggregated number of observed and total lineages
Kt_mean$TimePoints <- timePoints

## Melt Kt_mean for plotting
Kt_mean_melt <- melt(Kt_mean, id.vars = c("TimePoints", "NbObservedLin", "aggregNbHiddenLin", "aggregNbTotalLin", "Cred0.025", "Cred0.5", "Cred0.975"), variable.name = "NbHiddenLin", value.name = "ProbabilityDensity")
Kt_mean_melt$NbHiddenLin <- as.integer(Kt_mean_melt$NbHiddenLin)-1
# Kt_mean_melt <- Kt_mean_melt[Kt_mean_melt$ProbabilityDensity!=0,]                     # Remove 0-probability rows

## Get the distribution of the total number of lineages
Kt_mean_melt$NbTotalLin <- Kt_mean_melt$NbObservedLin + Kt_mean_melt$NbHiddenLin
Kt_mean_melt <- Kt_mean_melt[Kt_mean_melt$NbTotalLin < max(Kt_mean$Cred0.975)*1.1,]   # Remove lowest probability rows

## Plot densities and LTTs
ggplot(Kt_mean_melt, aes(x=TimePoints, y=NbHiddenLin, z = ProbabilityDensity)) + 
  geom_raster(aes(fill = ProbabilityDensity), interpolate = TRUE) +
  #geom_contour(colour = "white", binwidth = 0.002) +
  annotate(geom="line", x=Kt_mean$TimePoints, y=Kt_mean$aggregNbHiddenLin) +
  scale_fill_distiller(palette = "Spectral", direction = -1) +
  ggtitle("Probability density of the number of hidden lineages through time")

cols    <- c( "c1" = "dodgerblue3", "c2" = "gray25", "c3" = "forestgreen", "c95%" = "darkolivegreen4" )

ggplot(Kt_mean_melt, aes(x=TimePoints, y=NbTotalLin, z = ProbabilityDensity)) + 
  #geom_raster(aes(fill = ProbabilityDensity), interpolate = TRUE) +
  annotate(geom="raster", x=Kt_mean_melt$TimePoints, y=Kt_mean_melt$NbHiddenLin, interpolate = TRUE,
           fill = scales::colour_ramp(c(colorRampPalette(c("transparent", "dodgerblue2", "dodgerblue3", "dodgerblue4", "black"))(N)))(Kt_mean_melt$ProbabilityDensity)) +
  annotate(geom="raster", x=Kt_mean_melt$TimePoints, y=Kt_mean_melt$NbTotalLin, interpolate = TRUE,
           fill = scales::colour_ramp(colorRampPalette(c("transparent", "green4", "forestgreen", "black"))(N))(Kt_mean_melt$ProbabilityDensity)) +
  geom_line(aes(y=aggregNbHiddenLin, color="c1"), size=0.7) +
  geom_line(aes(y=NbObservedLin, color="c2"), size=0.7) + 
  geom_line(aes(y=aggregNbTotalLin, color="c3"), size=0.7) + 
  geom_line(aes(y=Cred0.025, color="c95%"), alpha=0) +
  annotate(geom="line", x=Kt_mean$TimePoints, y=Kt_mean$Cred0.025, color=cols["c95%"], linetype="twodash", size=0.5) +
  annotate(geom="line", x=Kt_mean$TimePoints, y=Kt_mean$Cred0.975, color=cols["c95%"], linetype="twodash", size=0.5) +
  scale_color_manual(name = "Lineages", breaks = c("c1", "c2", "c3", "c95%"), 
                     values = cols, labels = c("Hidden", "Observed", "Total", "95% credence interval")) +
  scale_x_continuous(name = "Time", expand = c(0.01,0.01)) +   
  scale_y_continuous(name = "Number of lineages", expand = c(0.01,0.01)) + 
  theme(panel.background=element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Probability density of the total number of lineages through time")
