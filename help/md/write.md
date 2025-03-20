## name
write
## title
Write RevObject to file
## description
This function write values in a RevObject to a file specified by the user.
## details
This function accepts multiple RevObjects in the first arguments to be written to a file.
After this, users can specify the filename with a string which can include the directory path
to where the file should be made. Users can also specify whether to append or overwrite the file
using a boolean operator (default is false). Lastly, a seperator can be specified using a string
for specifying how to separate values in the RevObject (default is "").
## authors
The RevBayes Development Core Team
## see_also
writeDelimitedCharacterData
writeFasta
## example

# define RevObject to write 
    x <- matrix([[1, 1],[1, 1]])
# write to CSV file
    write(x, "/path/to/file/filenmae.csv", false, ",")
## references
