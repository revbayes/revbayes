## name
fnBoundaryCrosserRates
## title
Boundary-crosser speciation and extinction rates
## description
Foote's boundary-crosser rates for each interval, estimated from the taxa that cross an interval boundary.
## details
Taxa confined to a single interval are dropped, which removes much of the bias incomplete sampling introduces but discards data: the method improves as sampling improves, and does worst when sampling is poor. An interval with no taxa ranging through it has no estimate and returns NaN.

Follows Foote (2000) equations 22 and 23, in the form given by Warnock et al. (2020) equation 8.

The record enters as the taxa read with `readTaxonData`. A taxon's first and last appearance are the oldest and youngest intervals it was sampled in, which is how these methods read a binned record.

The value is a two by l matrix, speciation in the first row and extinction in the second, since both rates come from one pass over the same counts.

Intervals follow the range processes: `timeline` holds the rate shift times youngest first and `present` the minimum age, so `timeline = v(10, 20, 30)` defines the intervals [0,10) [10,20) [20,30) and an unbounded oldest interval. An interval with no information returns NaN rather than an error, and the unbounded oldest interval never has a rate.

The oldest interval is unbounded unless `max_age` bounds it. These estimators divide a per-interval proportion by the interval's duration, so an unbounded interval has no rate. Occurrences older than `max_age` fall outside every interval and are not counted, which the function warns about once.

An occurrence whose reported bin straddles a boundary is counted at the interval holding the bin's midpoint, which is what the comparison literature does. `ambiguous="overlap"` counts it in every interval the bin touches and `ambiguous="exclude"` drops it.
## authors
June Walker
## see_also
fnPerCapitaRates
fnThreeTimerRates
## example
taxa <- readTaxonData("fossils.tsv")
timeline <- v(10, 20, 30)

rates := fnBoundaryCrosserRates(taxa, timeline)
speciation := rates[1]
extinction := rates[2]

# The estimates are fixed by the record, so log them once with mnFile rather than
# monitoring them.
monitors.append( mnFile(filename="output/rates.log", rates, printgen=1) )
## references
	- citation: Origination and extinction components of taxonomic diversity: general problems. Foote, Mike. 2000. Paleobiology, 26:74-102.
	  doi: https://doi.org/10.1017/S0094837300026890
	- citation: Assessing the impact of incomplete species sampling on estimates of speciation and extinction rates. Warnock, Rachel C. M., Heath, Tracy A. and Stadler, Tanja. 2020. Paleobiology, 46:137-157.
	  doi: https://doi.org/10.1017/pab.2020.12
