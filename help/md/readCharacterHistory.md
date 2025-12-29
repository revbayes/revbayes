## name
readCharacterHistory
## title
Function to read in character histories
## description
Reads character histories (in `simmap` format) from a file and saves them in one vector object. The character states always should be 0-indexed.
## details
A file name must be provided as argument.

`char_hist = readCharacterHistory(...)` returns a `CharacterHistory[]` object.
`char_hist = readCharacterHistory(...)[i]` returns a `CharacterHistory` object.


## authors
Priscilla Lau
## see_also
## example
    # read a character history
    char_hist = readCharacterHistory("output/simmap.tree")
    
## references
