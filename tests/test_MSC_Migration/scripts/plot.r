library(ggplot2)
library(dplyr)
library(tidyverse)

Ne <- c()
pow_list <- seq(1,5,0.2)

for (i in 1:length(pow_list)) {
    Ne[i] <- 0.0001 * (10^(pow_list[i]))
}

MSC <- c(-24.009, -13.858, -7.794, -4.307, -2.447, -1.614, -1.428, -1.650, -2.130, -2.773, -3.519, -4.329, -5.180, -6.057, -6.951, -7.854, -8.764, -9.678, -10.595, -11.513, -12.432)

MSC_M <- c(-23.965, -13.858, -7.794, -4.307, -2.447, -1.614, -1.428, -1.650, -2.130, -2.773, -3.519, -4.329, -5.180, -6.057, -6.951, -7.854, -8.764, -9.678, -10.595, -11.513, -12.432)

MSC_InverseGamma <- c(-9.956, -7.772, -5.844, -4.240, -3.015, -2.197, -1.777, -1.711, -1.935, -2.380, -2.984, -3.698, -4.485, -5.320, -6.187, -7.073, -7.972, -8.879, -9.791, -10.707, -11.624)

df <- data.frame(Ne, MSC, MSC_M, MSC_InverseGamma)

# Pivot data from wide to long
df <- df %>% 
  pivot_longer(cols = names(df)[2:4], values_to = "lnProbability", names_to = "model") 

df

# plot both likelihoods in one figure
ggplot(df, aes(x = Ne, y = lnProbability, col=model, label=round(lnProbability,1))) +
    geom_point(size=4) +
    geom_line() +
    scale_x_continuous(trans = 'log10') +
    geom_text(size=3, vjust = -2,angle = 45)

ggsave("tests/test_MSC_Migration/plots/lnProbability.png")

# plot only the msc-m
df_mscm <- df %>% filter(model == "MSC_M")
ggplot(df_mscm , aes(x = Ne, y = lnProbability, label=round(lnProbability,1))) +
    geom_point(size=4) +
    geom_line() +
    scale_x_continuous(trans = 'log10') +
    geom_text(size=3, vjust = -2,angle = 45)
ggsave("tests/test_MSC_Migration/plots/lnProbability_mscm.png")
