library(ggplot2)
library(dplyr)
library(tidyverse)

Ne <- c()
pow_list <- seq(1,5,0.2)

for (i in 1:length(pow_list)) {
    Ne[i] <- 0.0001 * (10^(pow_list[i]))
}

MSC <- c(-3822.248, -2406.986, -1514.357, -951.485, -596.677, -373.148, -232.451, -144.017, -88.559, -53.907, -32.383, -19.142, -11.128, -6.411, -3.775, -2.451, -1.956, -1.983, -2.341, -2.906, -3.603)

MSC_M <- c(-23.058, -22.557, -22.771, -23.146, -25.695, -24.636, -28.071, -28.975, -21.751, -21.593, -14.570, -8.422, -4.875, -2.978, -2.120, -1.919, -2.132, -2.606, -3.245, -3.988, -4.797)

df <- data.frame(Ne, MSC, MSC_M)

df

# Pivot data from wide to long
df <- df %>% 
  pivot_longer(cols = names(df)[2:3], values_to = "lnProbability", names_to = "model") 

# plot both likelihoods in one figure
ggplot(df, aes(x = Ne, y = lnProbability, col=model, label=round(lnProbability,1))) +
    geom_point(size=4) +
    geom_line() +
    scale_x_continuous(trans = 'log10') +
    geom_text(size=3, vjust = -2,angle = 45)

ggsave("02_likelihood_msc_vs_mscm/plots/lnProbability.png")

# plot only the msc-m
df_mscm <- df %>% filter(model == "MSC_M")
ggplot(df_mscm , aes(x = Ne, y = lnProbability, label=round(lnProbability,1))) +
    geom_point(size=4) +
    geom_line() +
    scale_x_continuous(trans = 'log10') +
    geom_text(size=3, vjust = -2,angle = 45)
ggsave("02_likelihood_msc_vs_mscm/plots/lnProbability_mscm.png")
