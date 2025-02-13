## name
readFossilCountsFile
## title
Reading fossil counts vectors and matrices
## description
This function is used to read fossil counts files for use with fossilized birth-death range models
## details
readFossilCountsFile takes on five parameters: file (the filename of the fossil counts matrix), header (whether to skip the first line, false my default), separator/delimiter (the character used to separate cells in the file, whitespace by default), rownames (whether to skip the first column, false by default), and taxa (taxa names for the fossil counts in the file). If the file contains a vector (i.e. only one line, considering rownames and header arguments), it will be interpreted as fossil counts for each interval. If it contains a matrix (i.e. multiple lines), it will be interpreted as fossil counts for each species (rows) for each interval. Note that the function expects taxon names (to match with the names contained in taxa), so that in the case where the file contains only a vector, the first column will be ignored (and therefore the file should contain one more column than the number of intervals).
## authors
Bruno do Rosario Petrucci, David Cerny
## see_also
readDataDelimitedFile
## example
## references
