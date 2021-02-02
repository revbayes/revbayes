## name
getOption
## title
Get a global RevBayes option
## description
Get a global option for RevBayes.
## details
Options are used to personalize RevBayes and are stored on the local machine. Currently this is rather experimental.
## authors
Sebastian Hoehna
## see_also
setOption
## example
	# compute the absolute value of a real number
	getOption("linewidth")
	
	# let us set the linewidth to a new value
	setOption("linewidth", 200)
	
	# now let's check what the value is
	getOption("linewidth")
	
## references
