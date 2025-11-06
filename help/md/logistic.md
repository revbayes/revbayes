## name
logistic
## title
The logistic function
## description
Compute the logistic function        
## details
The function is defined as

        logistic(x) = 1/(1 + exp(-x))

                    = exp(x)/(1 + exp(x))

This function takes a real number to a probability.
It is the inverse of the logit function.        
## authors
## see_also
logit        
## example
x ~ dnNormal(0,1)
p := logistic(x)                
## references
