## name
mvUpDownSlide
## title
Up-Down Sliding Proposal for several parameters jointly
## description
mvUpDownSlide move applies a sliding adjustment to multiple parameters simultaneously, where one parameter is increased while the other is decreased by the same amount
## details
This proposal randomly slides all a set of parameter up while the other set of parameters is slided down by the same value.
This should hopefully improve mixing in many cases.The actual sliding factor is computed by delta = lambda * ( u - 0.5 )
where u ~ Uniform(0,1)
## authors
## see_also
mvSlide
## example
p1 ~ dnUniform(0,1)
p2 ~ dnUniform(0,1)
moves.append(mvUpDownSlide(p1, p2, delta=0.05, weight=1))
## references
- citation: Yang, Ziheng, Molecular Evolution: A Statistical Approach (Oxford, 2014; online edn, Oxford Academic, 21 Aug. 2014), https://doi.org/10.1093/acprof:oso/9780199602605.001.0001, accessed 16 Apr. 2025.