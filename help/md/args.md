## name
args
## title
Vector of command-line arguments
## description
Vector holding command-line arguments supplied using the `--args` or `--cmd` flags. If no command-line arguments are supplied, `args` is initialized to an empty vector. Like a regular RevBayes vector, `args` uses 1-based indexing (i.e., the first element is accessed using `args[1]`), but unlike a regular vector, it can hold elements of different types.
## details
## authors
Sebastian Hoehna
## see_also
## example
	# args is initialized to an empty vector if no command-line arguments are supplied.
    # Assume RevBayes was called as follows: ./rb   
	args                          # returns [ ]
    args.size()                   # returns 0
    
    # However, args is not a reserved keyword. We can remove it from the workspace or overwrite it.
    clear(args)
    args <- readTrees("primates.tree")[1]  # args is now a tree
    
    # Assume RevBayes was called as follows: ./rb --args 1 5
    args[1] + args[2]             # returns 6
    
    # Unlike regular vectors, args can hold values of different types.
    # Assume RevBayes was called as follows: ./rb --args 1 5 2 "Hello " "world"
    print(args[4] + args[5])      # prints "Hello world"
    type(args[4])                 # returns String
    
    (args[1] + args[2])^args[3]   # returns 36
    type(args[1])                 # returns Natural
## references
