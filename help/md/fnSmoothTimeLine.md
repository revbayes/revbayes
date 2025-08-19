## name
fnSmoothTimeLine
## title
Create a smooth timeline
## description
Function to create a smooth timeline where all values after a maximum time are constant, i.e., equal to the previous interval, to avoid crazy looking plots from the prior.
## details
Thus function takes a vector of values and a matching vector of times and a maximum time. Then, it constructs a smooth timeline by using all values before the maximum, and replacing all values after the maximum with the last value before the maximum. Thus, the timeline is smooth after the maximum.
## authors
Sebastian Hoehna
## see_also

## example

## references
