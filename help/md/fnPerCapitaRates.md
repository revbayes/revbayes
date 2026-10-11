## name
fnPerCapitaRates
## title
Per-taxon speciation and extinction rates
## description
Foote's per-taxon rates for each interval: the proportion of taxa whose first or last appearance falls in that interval.
## details
This estimator does not account for incomplete sampling, and every taxon contributes, including those confined to a single interval.

Follows Foote (2000) equations 12 and 13, in the form given by Warnock et al. (2020) equation 7.

A taxon's first and last appearance are the oldest and youngest intervals it was sampled in. The value is a two by l matrix, speciation in the first row and extinction in the second.

`timeline` holds the rate shift times youngest first and `present` the minimum age, so `timeline = v(10, 20, 30)` defines the intervals [0,10) [10,20) [20,30) and an unbounded oldest interval. An interval whose rate is undefined returns NaN.

The oldest interval has no rate unless `max_age` bounds it. Occurrences older than `max_age` are not counted, and a `max_age` no older than the oldest rate shift leaves that interval NaN. Both cases warn once.

An occurrence whose reported bin straddles a boundary is counted at the interval holding the bin's midpoint. `ambiguous="overlap"` counts it in every interval the bin touches and `ambiguous="exclude"` drops it.
## authors
June Walker
## see_also
fnBoundaryCrosserRates
fnThreeTimerRates
## example
taxa <- readTaxonData("fossils.tsv")
timeline <- v(10, 20, 30)

rates := fnPerCapitaRates(taxa, timeline)
speciation := rates[1]
extinction := rates[2]

# mnModel skips a node that depends only on constants, so log it with mnFile
monitors.append( mnFile(filename="output/rates.log", rates, printgen=1) )
## references
	- citation: Origination and extinction components of taxonomic diversity: general problems. Foote, Mike. 2000. Paleobiology, 26:74-102.
	  doi: https://doi.org/10.1017/S0094837300026890
	- citation: Assessing the impact of incomplete species sampling on estimates of speciation and extinction rates. Warnock, Rachel C. M., Heath, Tracy A. and Stadler, Tanja. 2020. Paleobiology, 46:137-157.
	  doi: https://doi.org/10.1017/pab.2020.12
