## name
time
## title
Get the time information
## description
Get the current system time.

## details

"year" reports the current year (e.g. 2000).

"day" returns the index of the day in the year (e.g. Jan 1 = 1; Feb 1 = 32).

"(milli/nano)seconds" returns the number of (milli/nano)seconds that have elapsed since midnight.

"fromBeginning", the default, returns the number of milliseconds that have elapsed since 1400-Jan-01 00:00:00, the earliest representable date in the boost library's implementation of the Gregorian date system.

## authors
Sebastian Hoehna
## see_also
## example
	time()

	# Wait a little bit
	sum = 0
	for (i in 1:10000) sum += i
	# Now print the time again
	time()

## references
