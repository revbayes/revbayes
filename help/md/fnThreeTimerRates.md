## name
fnThreeTimerRates
## title
Three-timer speciation and extinction rates
## description
Alroy's three-timer rates for each interval, with the counts corrected by the sampling probability of a neighbouring interval.
## details
A taxon counts toward an interval only if it was sampled in that interval's neighbours as well, so the oldest and youngest intervals return NaN. An interval's sampling probability is its three-timer count over its three-timer plus part-timer count. Part-timers are taxa sampled on both sides of an interval but not within it. Speciation is corrected by the sampling probability of the next older interval and extinction by that of the next younger, so a rate is NaN where that neighbour has neither three-timers nor part-timers.

`extant=TRUE` reads the record's extant status, so a taxon reported as alive ranges to the present and crosses the younger edge of every interval it occupies, the youngest included. The default is FALSE, because a taxon whose status was never set reads as alive. `readTaxonData` takes the status from a status column, or from whether the taxon was sampled at the present. Warnock et al. (2020) assume complete knowledge of extant diversity for these methods.

The method assumes intervals of equal length, and warns once if they are not. A rate can come out negative, where the sampling correction exceeds the counts it was applied to, and is returned unchanged.

Follows the corrected formulas of Alroy (2008).

A taxon's first and last appearance are the oldest and youngest intervals it was sampled in. The value is a two by l matrix, speciation in the first row and extinction in the second.

`timeline` holds the rate shift times youngest first and `present` the minimum age, so `timeline = v(10, 20, 30)` defines the intervals [0,10) [10,20) [20,30) and an unbounded oldest interval. An interval whose rate is undefined returns NaN.

The oldest interval has no rate unless `max_age` bounds it. Occurrences older than `max_age` are not counted, and a `max_age` no older than the oldest rate shift leaves that interval NaN. Both cases warn once.

By default an occurrence whose reported bin straddles a boundary is dropped, keeping only occurrences resolved to a single interval. `ambiguous="midpoint"` counts it at the interval holding the bin's midpoint, `ambiguous="majority"` at the interval the bin overlaps most (the older one on a tie), and `ambiguous="all"` in every interval the bin touches. To resample ages within their bins, draw them in the script and call the function once per draw.
## authors
June Walker
## see_also
fnPerTaxonRates
fnBoundaryCrosserRates
## example
taxa <- readTaxonData("fossils.tsv")
timeline <- v(10, 20, 30)

rates := fnThreeTimerRates(taxa, timeline)
speciation := rates[1]
extinction := rates[2]

# mnModel skips a node that depends only on constants, so log it with mnFile
monitors.append( mnFile(filename="output/rates.log", rates, printgen=1) )
## references
	- citation: Dynamics of origination and extinction in the marine fossil record. Alroy, John. 2008. Proceedings of the National Academy of Sciences, 105:11536-11542.
	  doi: https://doi.org/10.1073/pnas.0802597105
	- citation: Assessing the impact of incomplete species sampling on estimates of speciation and extinction rates. Warnock, Rachel C. M., Heath, Tracy A. and Stadler, Tanja. 2020. Paleobiology, 46:137-157.
	  doi: https://doi.org/10.1017/pab.2020.12
