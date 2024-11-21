## name
Simplex
## title
Simplex
## description
A simplex is a vector of elements that sum to 1.
## details
## authors
## see_also

Moves that operate on Simplexes:
- mvBetaSimplex
- mvDirichletSimplex
- mvDPPValueBetaSimplex
- mvElementSwapSimplex

## example
```rb
x <- simplex([2, 2, 6])
x # = [ 0.2, 0.2, 0.6]
sum(x) # 1, by definition
```
## references
