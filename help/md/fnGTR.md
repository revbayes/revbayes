## name
fnGTR
## title
The General Time-Reversible rate matrix
## description
The GTR rate matrix.
## details
The general time-reversible rate matrix:
  Q(i,j) = exchangeRates(i,j) * baseFrequencies[j]
The exchangeRates are symmetric.
## authors
## see_also
## example
        er ~ dnDirichlet( v(1,1,1,1,1,1) )
        pi ~ dnDirichlet( v(1,1,1,1) )
        Q := fnGTR(er,pi)
## references
- citation: Tavare, S. Some Probabilistic and Statistical Problems in the Analysis of DNA Sequences.  Lectures on Mathematics in the Life Sciences (1986). 17: 57-86
  url: http://www.damtp.cam.ac.uk/user/st321/CV_&_Publications_files/STpapers-pdf/T86.pdf
  doi: null
