## name
fnThreeTimerRates
## title
Three-timer speciation and extinction rates
## description
Alroy's three-timer rates for each interval, with the counts corrected by a sampling probability pooled over the record.
## details
A taxon counts toward an interval only if it was sampled in that interval's neighbours as well, so the oldest and youngest intervals have no estimate and return NaN. The sampling probability is the three-timer count over the three-timer plus part-timer count, summed across every interval.

The method assumes intervals of equal length, and warns once if they are not.
Warnock's fbdR, the implementation behind Warnock et al. (2020), floors a negative three-timer rate at zero. This function returns it as it comes, since a negative estimate says the correction has outrun the counts and that is worth seeing.

Follows Alroy (2008), in the form given by Warnock et al. (2020) equations 9 and 10.

The record enters as the taxa read with `readTaxonData`. A taxon's first and last appearance are the oldest and youngest intervals it was sampled in, and the correction is estimated from part-timers, taxa sampled on both sides of an interval but not within it.

The value is a two by l matrix, speciation in the first row and extinction in the second, since both rates come from one pass over the same counts.

Intervals follow the range processes: `timeline` holds the rate shift times youngest first and `present` the minimum age, so `timeline = v(10, 20, 30)` defines the intervals [0,10) [10,20) [20,30) and an unbounded oldest interval. An interval with no information returns NaN rather than an error, and the unbounded oldest interval never has a rate.

The oldest interval is unbounded unless `max_age` bounds it. These estimators divide a per-interval proportion by the interval's duration, so an unbounded interval has no rate. Occurrences older than `max_age` fall outside every interval and are not counted, which the function warns about once.

An occurrence whose reported bin straddles a boundary is counted at the interval holding the bin's midpoint, which is what the comparison literature does. `ambiguous="overlap"` counts it in every interval the bin touches and `ambiguous="exclude"` drops it.
## authors
June Walker
## see_also
fnPerCapitaRates
fnBoundaryCrosserRates
## example
taxa <- readTaxonData("fossils.tsv")
timeline <- v(10, 20, 30)

rates := fnThreeTimerRates(taxa, timeline)
speciation := rates[1]
extinction := rates[2]

# The estimates are fixed by the record, so log them once with mnFile rather than
# monitoring them.
monitors.append( mnFile(filename="output/rates.log", rates, printgen=1) )
## references
	- citation: Dynamics of origination and extinction in the marine fossil record. Alroy, John. 2008. Proceedings of the National Academy of Sciences, 105:11536-11542.
	  doi: https://doi.org/10.1073/pnas.0802597105
	- citation: Assessing the impact of incomplete species sampling on estimates of speciation and extinction rates. Warnock, Rachel C. M., Heath, Tracy A. and Stadler, Tanja. 2020. Paleobiology, 46:137-157.
	  doi: https://doi.org/10.1017/pab.2020.12
