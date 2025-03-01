## name
readVCF
## title
Read VCF
## description
Read VCF file into RevBayes
## details
readVCF reads in a file that is in Variant Call Format (VCF), accepting two
 arguments. The first argument specifies the relative or absolute 
file path to desired VCF file. The second specifies type of data
to be constructed (default binary). This function
only allows for 0, 1, and . characters in the VCF file.
## authors
## see_also
## example
    x <- readVCF("path/to/VCF/file", "binary")
## references
Danecek, P., Auton, A., Abecasis, G., Albers, C. A., Banks, E., DePristo, M. A., ... & 1000 Genomes Project Analysis Group. (2011). The variant call format and VCFtools. Bioinformatics, 27(15), 2156-2158.

