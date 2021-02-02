## name
ifelse
## title
If-else statement as a function
## description
If the expression is true, then the function returns the first value, otherwise the second value.
## details
The ifelse function is important when the value of a variable should deterministically change during an analysis depending on other variables. Standard if-else statements are not dynamically re-evaluated.
## authors
Sebastian Hoehna
## see_also
## example
	a <- 1
	b := ifelse( a == 1, 10, -10 )
	b
	
	a <- 2
	b
	
## references
