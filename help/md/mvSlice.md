## name
mvSlice
## title
Propose a slice move
## description
`mvSlice` proposes a new value for a variable based on the current shape of its likelihood function.

## details
A slice proposal uses the shape of the current likelihood distribution of a variable to propose a new value.
First, a likelihood value is drawn uniformly in order to define a horizontal 'slice' through the likelihood distribution.
Then, a new value is drawn uniformly from those values that lie within this slice.

This allows parameter space to be traversed more efficiently than a random walk.
A practical outcome of the implementation is that small moves are proposed in certain parts of parameter space, and large moves in other parts of the space, as appropriate.

A detailed explanation with figures is provided in Neal (2003).

## authors
## see_also
`mvSlide` and `mvScale` are possible alternatives where a fixed move size is desired.
## example
## references
	- citation: Slice sampling. Neal (2003). Ann. Statist. 31(3): 705-767    
	  doi: 10.1214/aos/1056562461
	  url: https://projecteuclid.org/journals/annals-of-statistics/volume-31/issue-3/Slice-sampling/10.1214/aos/1056562461.full
