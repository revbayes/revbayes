## name
fnJC
## title
fnJC
## description
Jukes Cantor rate matrix
## details
Rate matrix with n states in which each state has an equal probability to change 
to any other state. The rate of transition from one state to another is equal to 
n / n-1.
## authors
Sebastian Hoehna 
## see_also
fnGTR
fnF81
## example
    # Rate Matrix for DNA
    q := fnJC(4)
    # Rate Matrix for Amino Acids
    q := fnJC(20)
    # Binary Character Matrix
    q := fnJC(2)
## references
- citation: Jukes TH, Cantor CR (1969). Evolution of Protein Molecules. New York: Academic Press. pp. 21–132.
