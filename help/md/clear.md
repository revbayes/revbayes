## name
clear
## title
Clear the current workspace
## description
Clear (e.g., remove) variables and functions from the workspace.
## details
The clear function removes either a given variable or all variables from the workspace. Clearing the workspace is very useful between analysis if you do not want to have old connections between variables hanging around.
## authors
Sebastian Hoehna
## see_also
exists
## example
	ls()   # check what is in the workspace
	a <- 1
	b := exp(a)
	ls()   # check what is in the workspace
	clear()
	ls()   # check what is in the workspace
	a <- 1
	b := exp(a)
	ls()   # check what is in the workspace
	clear( b )
	ls()   # check what is in the workspace
	
## references
