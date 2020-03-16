## name
sort
## title
Sort function
## description
Function for sorting the members of a vector in either ascending or descending order.
## details
The vector to be sorted can be of any numeric type. Ascending or descending is specified via the `ascending` argument
## authors
## see_also
## example
    nums = v(1,3,5,7,2,4,6,8)
    sort(nums)
    # this will result in 1,2,3,4,5,6,7,8
    sort(nums, ascending = FALSE) 
    # this will result in 8,7,6,5,4,3,2,1
## references
