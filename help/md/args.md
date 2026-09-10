## name
args
## title
Vector of command-line arguments
## description
When RevBayes is called from the command line, any tokens following the name
of a script file (or an `-e` expression) are treated as arguments to the script
(or to the expression) and stored in the `args` vector.
## details
If no command-line arguments are supplied, `args` is initialized to an empty
vector. Like regular RevBayes vectors, `args` uses 1-based indexing (i.e., the
first element is accessed using `args[1]`), but unlike a regular vector, it can
hold elements of different types.
## authors
Sebastian Hoehna
Ben Redelings
## see_also
## example
    # args is initialized to an empty vector if no command-line arguments are supplied.
    # Assume RevBayes was called as follows: ./rb   
    args                          # returns [ ]
    args.size()                   # returns 0
    
    # However, args is not a reserved keyword. We can remove it from the workspace or overwrite it:
    clear(args)
    args <- readTrees("primates.tree")[1]  # args is now a tree
    
    # To illustrate how args works, run RevBayes in interactive mode with an empty expression:
    # ./rb -i -e "" 1 5
    args[1] + args[2]             # returns 6
    
    # Unlike regular vectors, args can hold values of different types. Call RevBayes as follows:
    # ./rb -i -e "" 1 5 2 "Hello " "world"
    print(args[4] + args[5])      # prints "Hello world"
    type(args[4])                 # returns String
    
    (args[1] + args[2])^args[3]   # returns 36
    type(args[1])                 # returns Natural
## references
