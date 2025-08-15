## name
setOption
## title
Set a global RevBayes option
## description
Set a global option for RevBayes.
## details
Runtime options are used to personalize RevBayes and are stored on the local machine.
The currently available keys and their associated values are as follows:

    lineWidth=<integer>          Screen width when printing (in characters).

        DEFAULT: 160

    outputPrecision=<integer>    How many significant digits to print for the values of model graph nodes.

        DEFAULT: 7

    printNodeIndex=<TRUE,FALSE>  Print the node indices of a tree as annotations?

        DEFAULT: TRUE

    useScaling=<TRUE,FALSE>      Should the partial likelihoods in continuous-time Markov chain (CTMC) models be scaled
                                 to avoid underflow?
        DEFAULT: TRUE

    scalingDensity=<integer>     If so, scale CTMC likelihoods every n-th node (min = 1).

        DEFAULT: 1

    tolerance=<numeric>          Tolerance for comparing doubles.

        DEFAULT: 10e-10

    debugMCMC=<0,1>              How much work to perform to check MCMC?

        0: MCMC run without checks.
        1: MCMC run with additional checks at extra CPU time cost.

        DEFAULT: 0

    logMCMC=<0,1,2,3,4>          How much logging to perform when checking MCMC? NOTE: This option serves for debugging and
                                 should be considered experimental. The exact meaning of individual values may be subject
                                 to frequent changes.

        0: No information on individual moves written out.
        1 or higher: Writes out the generation, within-generation position, and name for each move.
        2 or higher: Also writes out posterior, likelihood, prior, and Hastings ratios, and if log likelihood = -Inf or NaN,
                     writes out why this is the case.
        3 or higher: Writes out each changed probability density and the name of the corresponding model graph node.
        4: Writes out additional details about the mvSlice move (if present).

        DEFAULT: 0

    interactive=<TRUE,FALSE>     Read commands from STDIN and evaluate them until q() is received.

        DEFAULT: TRUE if no script or -e expr is given, otherwise FALSE.

    echo=<TRUE,FALSE>            Should commands be printed to the screen?

        DEFAULT: TRUE if interactive, otherwise FALSE.

    errorExit=<TRUE,FALSE>       Should we exit on the first error?

        DEFAULT: FALSE if interactive, otherwise otherwise TRUE.

## authors
Sebastian Hoehna
## see_also
getOption
## example
	# compute the absolute value of a real number
	getOption("linewidth")
	
	# let us set the linewidth to a new value
	setOption("linewidth", 200)
	
	# now let's check what the value is
	getOption("linewidth")
	
## references
