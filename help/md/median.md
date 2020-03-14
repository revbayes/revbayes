## name
median
## title
Median of a set of numbers
## description
Finds the median of a sorted vector of numbers.
## details
The vector is sorted when `median` is used finding the
number of the sorted values with an equal amount of numbers that
are greater than or less than that value. If the length of the vector is even, there will be no such value. In that case, the two are averaged automatically.
## authors
## see_also
`mean`
## example
    a = v(5,3,2,6,8)
    median(a)
    # 5 is the result
    b = v(1,1,2,3,5,8)
    median(b)
    # 2.5 is the result
## references
