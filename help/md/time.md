## name
time
## title
Get the time information
## description
Get the current system time.

## details

"year" reports the current year (e.g. 2000).

"day" returns the index of the day in the year (e.g. Jan 1 = 1; Feb 1 = 32).

"(milli)seconds" returns the number of (milli)seconds that have elapsed since midnight.

"fromBeginning", the default, returns the number of milliseconds that have elapsed since 1 January 1970, the start of the Unix epoch.

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
