## name
fnJC
## title
fnJC
## description
Jukes cantor transition matrix
## details
Transition matrix with n amount of states in which each state has an equal
probability of to change to any other state. The only parameter in this model
is mu with the rate of change from any given state being mu / number of states.
## authors
Jukes & Cantor 
## see_also
tbd ~ other mutation models?
## example
    # Rate Matrix for DNA
    q := fnJC(4)
    # Rate Matrix for Amino Acids
    q := fnJC(20)
    # Binary Character Matrix
    q := fnJC(2)
## references
Jukes TH, Cantor CR (1969). Evolution of Protein Molecules. New York: Academic Press. pp. 21–132.
