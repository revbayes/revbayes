## name
exp
## title
Exponential of a number
## description
Maps the value of a number x to e^x, where e is the number such that `ln(e) = 1`.
## details
## authors
## see_also
## example
	# checking that ln(e) = 1
	x <- exp(1)
	ln_of_x <- ln(x)
	if (ln_of_ex != 1) {
		print("Problem when computing an exponential value.")
	} else {
		print("Correct computation of an exponential value.")
	}
## references
