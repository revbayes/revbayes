## name
CharacterHistory
## title
Discrete character history
## description
Stores a discrete character history, usually produced by reading in a file in `simmap` format, or retrieving from a `dnPhyloCTMCDASiteIID` distribution.
## details
Individual character histories can be accessed as follows:

- `char_hist[n]`: returns the `n`th character history.

## authors
Priscilla Lau
## see_also
## example
    # read character histories
    char_hist = readCharacterHistory("simmap.tree")
    
    # retrieve the first character history
    first_char_hist = char_hist[1]

## references
