setwd("Documents/Scolaire/ENS/Stage_M1/revbayes/tests/test_OBDP/output_inferKT_macro_test_cluster/")

devtools::load_all("../../../../RevGadgets")

Kt_mean <- rev.process.nbLineages( popSize_distribution_matrices_file="Kt_trace.txt", 
                                   trees_trace_file="../output_macro_test_Mt/mcmc_OBDP_macro_test.trees", 
                                   weight_trees_posterior=TRUE )

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
