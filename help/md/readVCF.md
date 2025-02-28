## name
readVCF
## title
Read VCF
## description
Read VCF file into RevBayes
## details
readVCF reads in a file that is in Variant Call Formati (VCF), accepting two
 arguments. The first argument specifies the relative or absolute 
file path to desired VCF file. The second specifies type of data
to be constructed (default binary). This function
only allows for 0, 1, and . characters in the VCF file.
## authors
## see_also
## example
    x <- readVCF("path/to/VCF/file", "binary")
## references
