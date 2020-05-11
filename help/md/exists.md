## name
exists
## title
Check whether a variable exists
## description
Determines whether the RevBayes workspace contains a variable named 'name'
## details
'exists' returns 'true' if the workspace contains a variable whose name matches the String 'name' and 'false' otherwise. One use of 'exists' is to add Move and Monitor objects conditional on the variable 'x' existing. The function 'ls' provides a summary for all variable names that 'exists' would evaluate as 'true'.
## authors
Michael Landis
## see_also
clear
## example
	## Correct usage: does "x" exist?
	x <- 1.0
	exists("x")
	
	## Incorrect usage: does "1.0" exist?
	exists(x)
	
## references
