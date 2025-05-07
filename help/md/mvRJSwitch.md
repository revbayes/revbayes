## name
mvRJSwitch
## title
A Reversible-jump swich move for parameter inclusion 
## description
The mvRJSwitch move implements a reversible-jump Markove Chain Monte Carlo (RJ-MCMC) proposal that allows a model parameter to swich between being active an inactive 
## details
The mvRJSwitch move enables reversible-jump MCMC (RJ-MCMC) transitions between two regimes. The move paired with dnRJMixture distribution which defines a prior that allows a parameter to either: take a fixed value of 0 (rj_null_value=0) with probability "1-p ", indicating that the parameter is inactive (excluded from the model), or be drawn from a continuous distribution (e.g., dnNormal, dnExponential) with probability p, indicating that the parameter is active (included in the model).
## authors
## see_also
## example
#set up priors for feature effects
rj_null_value <- 0.0          # fixed "off-value" for RJMCMC
rj_prob       <- 0.5          # prob. of RJMCMC taking "off-value"

#prior of "on-value" for RJMCMC
bound <- 2
rj_base_sym_dist = dnUniform(-bound, bound)
rj_sym_dist = dnRJMixture(rj_null_value, rj_base_sym_dist, p=rj_prob)
#initialize categorical feature effects, create moves, add monitor variables
for (i in 1:feature_CW.size()) {
    sigma_w[i].setValue(0)
    moves.append( mvScale(sigma_w[i], weight=2) )
    moves.append( mvSlide(sigma_w[i], weight=2) )
    moves.append( mvRJSwitch(sigma_w[i], weight=3) )
    use_sigma_w[i] := ifelse(sigma_w[i] == 0.0, 0, 1)
}
## references
- citation: Landis M.J., Quintero I., Muñoz M.M., Zapata F., Donoghue M.J. 2022. Phylogenetic inference of where species spread or split across barriers. Proceedings of the National Academy of Sciences. 119.
