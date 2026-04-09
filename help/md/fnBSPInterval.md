## name
fnBSPInterval
## title
Expand a vector of values for Bayesian skyline method input
## description
A Bayesian skyline approach, as suggested by Drummond et al. (2005) and
Billenstein & Höhna (2024), requires as an input a vector of values. These
values can be the same for several consecutive entries. For example, if the
total number of intervals is 4, and we have the two values, 0.2 and 0.5, where
the first value is used once and the second used three times, we could create
a Bayesian skyline vector as `[0.2, 0.5, 0.5, 0.5]`. To construct such
a vector, we use the function `fnBSPInterval`.
## details
The function `fnBSPInterval` takes in two arguments: `x`, the values; and `n`,
the number of times each value is replicated. Currently, the function assumes
that all values are positive real numbers, as `fnBSPInterval` is used for
population sizes in coalescent methods.
## authors
Sebastian Höhna
## see_also
dnCoalescent
## example
    a <- [0.1, 0.2, 0.3]
    b <- [2, 1, 4]
    fnBSPInterval( a, b )
## references
- citation: Billenstein RJ, Höhna S (2024). Comparison of Bayesian coalescent skyline plot models for inferring demographic histories. Molecular Biology and Evolution, 41(5):msae073.
  doi: 10.1093/molbev/msae073
  url: https://academic.oup.com/mbe/article/41/5/msae073/7648822
- citation: Drummond AJ, Rambaut A, Shapiro B, Pybus OG (2005). Bayesian coalescent inference of past population dynamics from molecular sequences. Molecular Biology and Evolution, 22(5):1185--1192.
  doi: 10.1093/molbev/msi103
  url: https://academic.oup.com/mbe/article/22/5/1185/1066885
