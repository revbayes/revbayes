## name
mvDirichletSimplex
## title
Dirichlet Simplex move
## description
A Dirichlet-simplex proposal randomly changes some values of a [Simplex](https://revbayes.github.io/documentation/Simplex.html)
(a vector whose elements sum to 1). The other values change too because of renormalization.
 
First, some random indices are drawn.
Then, the proposal draws a new simplex `u ~ Dirichlet(val[index] * alpha)`, where alpha is the tuning parameter.
The new value is set to `u`.
The simplex is then renormalized.

## details
## authors
## see_also
- mvBetaSimplex
- mvDPPValueBetaSimplex
- mvElementSwapSimplex

## example
Usage examples can be found at https://revbayes.github.io/tutorials/morph_tree/V2.html and https://revbayes.github.io/tutorials/morph_ase/ase_free.html
## references
