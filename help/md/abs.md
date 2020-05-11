## name
abs
## title
Absolute value of a number
## description
The 'abs' function returns the absolute value of a number.
## details
## authors
Sebastian Hoehna
## see_also
ceil
floor
round
## example
	# compute the absolute value of a real number
	number <- -3.0
	absoluteValueOfTheNumber <- abs(number)
	if (number + absoluteValueOfTheNumber != 0.0) {
	    print("Problem when computing an absolute value.")
	} else {
	    print("Correct computation of an absolute value.")
	}
	
## references
