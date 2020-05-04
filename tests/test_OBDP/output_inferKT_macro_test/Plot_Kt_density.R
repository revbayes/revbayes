setwd("~/Documents/Scolaire/ENS/Stage_M1/revbayes/tests/test_OBDP/output_inferKt_macro_test/")

devtools::load_all("../../../../RevGadgets")

Kt_mean <- rev.process.nbLineages( start_time_trace_file = "start_time_trace.txt", 
                                   popSize_distribution_matrices_file="Kt_trace.txt", 
                                   trees_trace_file="../output_macro_test_Mt/mcmc_OBDP_macro_test.trees", 
                                   weight_trees_posterior=FALSE )

p <- rev.plot.nbLineages( Kt_mean,
                          xlab="Time",
                          ylab="Number of lineages",
                          col.Hidden = "dodgerblue3",
                          col.Observed = "gray25",
                          col.Total = "forestgreen",
                          col.Hidden.interval = "dodgerblue2",
                          col.Total.interval = "darkolivegreen4",
                          palette.Hidden = c("transparent", "dodgerblue2", "dodgerblue3", "dodgerblue4", "black"),
                          palette.Total = c("transparent", "green4", "forestgreen", "black"),
                          line.size=0.7,
                          interval.line.size=0.5,
                          show.Hidden=TRUE,
                          show.Observed=TRUE,
                          show.Total=TRUE,
                          show.intervals=TRUE,
                          show.densities=TRUE,
                          show.expectations=TRUE,
                          use.interpolate=TRUE )
p

# Compare to the real number of lineages
truth_times <- t(read.table("../output_macro_test_Mt_truth/truth_popSize_times.csv", sep=";"))
truth_nbLin <- t(read.table("../output_macro_test_Mt_truth/truth_popSize_number.csv", sep=";"))
truth <- data.frame(times=-truth_times, nbTotalLin=truth_nbLin)

p + annotate(geom="line", x=truth$times, y=truth$nbTotalLin, color="red")

# Compare to the real number of hidden lineages
real_tree <- read.tree("../output_macro_test_Mt_truth/truth_tree.tree")
real_tree <- collapse.singles(real_tree)                                 # Remove sampled ancestors
real_tree_LTT <- data.frame(ltt.plot.coords(real_tree))

getNbObservedLin <- function(t){
  if (t < real_tree_LTT$time[1]) return (0)
  first_consecutive_time_point <- which(real_tree_LTT$time >= min(t, 0))[1]
  return (real_tree_LTT$N[first_consecutive_time_point])
}
truth$nbObservedLin <- sapply(truth$times, getNbObservedLin)
truth$nbHiddenLin <- truth$nbTotalLin - truth$nbObservedLin

p + annotate(geom="line", x=truth$times, y=truth$nbHiddenLin, color="red")
