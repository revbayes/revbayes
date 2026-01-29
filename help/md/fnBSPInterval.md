## name
fnBSPInterval
## title
Expand value vector for Bayesian skyline method input
## description
A Bayesian skyline approach, as suggested by Drummond et al. (2005; https://doi.org/10.1093/molbev/msi103) and Billenstein and Höhna (2024; https://doi.org/10.1093/molbev/msae073), requires as an input a vector of values. These values can be the same for several consecutive entries. For example, if the total number of intervals is 4, and we have the two values 0.2 and 0.5, where the first value is used once and the second used three times, we could create a Bayesian skyline vector as `[0.2,0.5,0.5,0.5]`. To construct such a vector, we use the function `fnBSPInterval`.
## details
The function `fnBSPInterval` takes in two arguments: `x` which are the values, and `n` which are how often each value is replicated. Currently, this function assumes that all values are positive real numbers, as the `fnBSPInterval` is used for population sizes in coalescent methods.
## authors
Sebastian Höhna
## see_also
dnCoalescent
## example
a <- [0.1,0.2,0.3]
b <- [2,1,4]
fnBSPInterval( a, b )
## references
