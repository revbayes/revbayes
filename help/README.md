# RevBayes Help Database

The Markdown files stored in the `md` directory are used to generate the database of help entries that populate the help system in RevBayes.

The build embeds these files in the `rb` executable unchanged, and `src/core/help/RbHelpDatabase.cpp` parses them the first time help is requested. Editing a file here and rebuilding is all that is needed to update the help.

The build also yields the `rb-help2yml` executable, which is used to generate an exhaustive YAML database of all information contained in the RevBayes help system, including the information contained in the Markdown database, as well as hardcoded elements such as function arguments and return types. A web portal for this database is presented on the [RevBayes documentation website](https://revbayes.github.io/documentation/).

The fields in each Markdown help entry file are:

## name
This field is the name of the RevBayes type to which the entry is mapped.
## title
A brief name or title for the entry.
## description
A longer description of the RevBayes type, what it does or what it is used for.
## details
A more detailed theoretical description if applicable.
## authors
A list of authors, one per line.
## see_also
A list of other relevant RevBayes types, one per line.
## example
A (optionally indented) literal block of Rev code that demonstrates the use of the type.
## references
A list of references as YAML entries with the tags `citation`, `doi` and `url`.
Here is an example:

	- citation: Bayes Estimation of Species Divergence Times and Ancestral Population
	    Sizes Using DNA Sequences From Multiple Loci. Bruce Rannala and Ziheng Yang. GENETICS
	    August 1, 2003 vol. 164 no. 4 1645-1656.
	  doi: null
	  url: http://www.genetics.org/content/164/4/1645.short
	- citation: 'Bayesian Inference of Species Trees from Multilocus Data. Heled and Drummond
	    Mol. Biol Evol. 27 (3): 570-580, 2010.'
	  doi: 'DOI: https://doi.org/10.1093/molbev/msp274'
	  url: https://academic.oup.com/mbe/article/27/3/570/999753/Bayesian-Inference-of-Species-Trees-from
