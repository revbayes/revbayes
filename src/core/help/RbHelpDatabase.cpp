/**
 * This file was generated automatically.
 * It is not intended to be human readable.
 * See help/README.md for details.
 */

#include "RbHelpDatabase.h"

using namespace std;

RevBayesCore::RbHelpDatabase::RbHelpDatabase()
{
	help_strings[string("AbstractHomologousDiscreteCharacterData")][string("name")] = string(R"(AbstractHomologousDiscreteCharacterData)");
	help_arrays[string("Bool")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("Bool")][string("description")] = string(R"(Bool variables can be either `true` or `false` (`TRUE` or `FALSE` also work).)");
	help_strings[string("Bool")][string("example")] = string(R"(a <- FALSE
if(!a)
    print("a is not true")
# this will print the statement in parentheses)");
	help_strings[string("Bool")][string("name")] = string(R"(Bool)");
	help_strings[string("Bool")][string("title")] = string(R"(Datatype for logical variables.)");
	help_arrays[string("BootstrapAnalysis")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("BootstrapAnalysis")][string("name")] = string(R"(BootstrapAnalysis)");
	help_strings[string("BranchLengthTree")][string("description")] = string(R"(The Tree datatype stores information to describe the shared ancestry of a taxon set. Information includes taxon labels, topology, nodecount, and branch lengths. Tree objects also possess several useful methods to traverse and manipulate the Tree's value.)");
	help_strings[string("BranchLengthTree")][string("name")] = string(R"(BranchLengthTree)");
	help_arrays[string("BranchLengthTree")][string("see_also")].push_back(string(R"(TimeTree)"));
	help_arrays[string("BranchLengthTree")][string("see_also")].push_back(string(R"(BranchLengthTree)"));
	help_strings[string("BranchLengthTree")][string("title")] = string(R"(Tree datatype)");
	help_strings[string("CharacterHistoryRateModifier")][string("name")] = string(R"(CharacterHistoryRateModifier)");
	help_strings[string("CladogeneticProbabilityMatrix")][string("name")] = string(R"(CladogeneticProbabilityMatrix)");
	help_strings[string("CladogeneticSpeciationRateMatrix")][string("name")] = string(R"(CladogeneticSpeciationRateMatrix)");
	help_strings[string("ContinuousCharacterData")][string("name")] = string(R"(ContinuousCharacterData)");
	help_strings[string("CorrespondenceAnalysis")][string("name")] = string(R"(CorrespondenceAnalysis)");
	help_strings[string("DistanceMatrix")][string("name")] = string(R"(DistanceMatrix)");
	help_arrays[string("HillClimber")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("HillClimber")][string("description")] = string(R"(The HillClimber analysis object keeps a model and the associated moves and monitors. The object is used to run Markov chain Monte Carlo (HillClimber) simulation on the model, using the provided moves, to obtain a sample of the posterior probability distribution. During the analysis, the monitors are responsible for sampling model parameters of interest.)");
	help_strings[string("HillClimber")][string("details")] = string(R"( The HillClimber analysis object produced by a call to this function keeps copies of the model and the associated moves and monitors. The HillClimber analysis object is used to run Markov chain Monte Carlo (HillClimber) simulation on the model, using the provided moves, to obtain a sample of the posterior probability distribution. During the analysis, the monitors are responsible for sampling model parameters of interest.)");
	help_strings[string("HillClimber")][string("example")] = string(R"(# Create a simple model (unclamped)
a ~ exponential(1)
mymodel = model(a)

# Create a move vector and a monitor vector
moves[1] = mvScale(a, lambda=1.0, weight=1.0)
monitors[1] = mnFile(a,"output/out.log")

# Create an HillClimber object
myHillClimberObject = HillClimber( mymodel, monitors, moves)

# Run a short analysis
myHillClimberObject.burnin( generations = 400, tuningInterval = 100)
myHillClimberObject.run( generations = 400)

# print the summary of the operators (now tuned)
myHillClimberObject.operatorSummary())");
	help_strings[string("HillClimber")][string("name")] = string(R"(HillClimber)");
	help_arrays[string("HillClimber")][string("see_also")].push_back(string(R"(SimulatedAnnealing)"));
	help_strings[string("HillClimber")][string("title")] = string(R"(Hill-Climber analysis object)");
	help_strings[string("Integer")][string("name")] = string(R"(Integer)");
	help_strings[string("MatrixReal")][string("name")] = string(R"(MatrixReal)");
	help_strings[string("MatrixRealPos")][string("name")] = string(R"(MatrixRealPos)");
	help_strings[string("MatrixRealSymmetric")][string("name")] = string(R"(MatrixRealSymmetric)");
	help_strings[string("Natural")][string("name")] = string(R"(Natural)");
	help_arrays[string("Probability")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("Probability")][string("description")] = string(R"(A Probability is a real value between 0.0 and 1.0)");
	help_strings[string("Probability")][string("example")] = string(R"(# Create a RealPos
x <- 12/13
type(x)

# Convert the RealPos to Probability
y := Probability(x)
type(y))");
	help_strings[string("Probability")][string("name")] = string(R"(Probability)");
	help_strings[string("RateGenerator")][string("name")] = string(R"(RateGenerator)");
	help_strings[string("Real")][string("description")] = string(R"(The real number data type can hold any real number value.
Not to be confused with integers which are whole numbers, or
`Natural` which are the counting numbers (e.g. 1,2,3,...).)");
	help_strings[string("Real")][string("example")] = string(R"(a = 1.1
b = 1.3
c = a + b
# c will be 2.4)");
	help_strings[string("Real")][string("name")] = string(R"(Real)");
	help_arrays[string("Real")][string("see_also")].push_back(string(R"(`RealPos`, `Integer`, `Natural`, `Probability`)"));
	help_strings[string("Real")][string("title")] = string(R"(Real number data type)");
	help_strings[string("RealPos")][string("name")] = string(R"(RealPos)");
	help_strings[string("RevObject")][string("name")] = string(R"(RevObject)");
	help_strings[string("Simplex")][string("name")] = string(R"(Simplex)");
	help_strings[string("String")][string("name")] = string(R"(String)");
	help_strings[string("TimeTree")][string("description")] = string(R"(The Tree datatype stores information to describe the shared ancestryof a taxon set. Information includes taxon labels, topology, nodecount, and branch lengths. Tree objects also possess several usefulmethods to traverse and manipulate the Tree's value.)");
	help_strings[string("TimeTree")][string("name")] = string(R"(TimeTree)");
	help_arrays[string("TimeTree")][string("see_also")].push_back(string(R"(TimeTree)"));
	help_arrays[string("TimeTree")][string("see_also")].push_back(string(R"(BranchLengthTree)"));
	help_strings[string("TimeTree")][string("title")] = string(R"(Tree datatype)");
	help_strings[string("Tree")][string("description")] = string(R"(The Tree datatype stores information to describe the shared ancestryof a taxon set. Information includes taxon labels, topology, nodecount, and branch lengths. Tree objects also possess several usefulmethods to traverse and manipulate the Tree's value.)");
	help_strings[string("Tree")][string("name")] = string(R"(Tree)");
	help_arrays[string("Tree")][string("see_also")].push_back(string(R"(TimeTree)"));
	help_arrays[string("Tree")][string("see_also")].push_back(string(R"(BranchLengthTree)"));
	help_strings[string("Tree")][string("title")] = string(R"(Tree datatype)");
	help_strings[string("VectorMonitors")][string("name")] = string(R"(VectorMonitors)");
	help_strings[string("VectorMoves")][string("name")] = string(R"(VectorMoves)");
	help_strings[string("[]")][string("name")] = string(R"([])");
	help_arrays[string("abs")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("abs")][string("description")] = string(R"(The 'abs' function returns the absolute value of a number.)");
	help_strings[string("abs")][string("example")] = string(R"(# compute the absolute value of a real number
number <- -3.0
absoluteValueOfTheNumber <- abs(number)
if (number + absoluteValueOfTheNumber != 0.0) {
    print("Problem when computing an absolute value.")
} else {
    print("Correct computation of an absolute value.")
})");
	help_strings[string("abs")][string("name")] = string(R"(abs)");
	help_arrays[string("abs")][string("see_also")].push_back(string(R"(ceil)"));
	help_arrays[string("abs")][string("see_also")].push_back(string(R"(floor)"));
	help_arrays[string("abs")][string("see_also")].push_back(string(R"(round)"));
	help_strings[string("abs")][string("title")] = string(R"(Absolute value of a number)");
	help_strings[string("ancestralStateTree")][string("name")] = string(R"(ancestralStateTree)");
	help_strings[string("annotateTree")][string("name")] = string(R"(annotateTree)");
	help_arrays[string("append")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("append")][string("description")] = string(R"('append' adds an element to a vector.)");
	help_strings[string("append")][string("details")] = string(R"('append' creates a new vector that is the original vector plus the extra element.)");
	help_strings[string("append")][string("example")] = string(R"(a <- 1:3
b <- 4
c := append(a,b))");
	help_strings[string("append")][string("name")] = string(R"(append)");
	help_arrays[string("append")][string("see_also")].push_back(string(R"(rep)"));
	help_strings[string("append")][string("title")] = string(R"(Append a value)");
	help_strings[string("beca")][string("name")] = string(R"(beca)");
	help_strings[string("branchScoreDistance")][string("name")] = string(R"(branchScoreDistance)");
	help_arrays[string("ceil")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("ceil")][string("description")] = string(R"(The 'ceil' function maps the value of a number to the smallest following integer.)");
	help_strings[string("ceil")][string("example")] = string(R"(# compute the ceiling of a real number
number <- 3.4
ceiled_number <- ceil(number)
if (ceiled_number != 4.0) {
    print("Problem when computing a ceiled value.")
} else {
    print("Correct computation of a ceiled value.")
})");
	help_strings[string("ceil")][string("name")] = string(R"(ceil)");
	help_arrays[string("ceil")][string("see_also")].push_back(string(R"(abs)"));
	help_arrays[string("ceil")][string("see_also")].push_back(string(R"(floor)"));
	help_arrays[string("ceil")][string("see_also")].push_back(string(R"(round)"));
	help_strings[string("ceil")][string("title")] = string(R"(Ceiling value of a number)");
	help_strings[string("characterMapTree")][string("name")] = string(R"(characterMapTree)");
	help_strings[string("checkNodeOrderConstraints")][string("name")] = string(R"(checkNodeOrderConstraints)");
	help_arrays[string("choose")][string("authors")].push_back(string(R"(Michael Landis)"));
	help_strings[string("choose")][string("description")] = string(R"(Rev function to calculate the binomial coefficients.)");
	help_strings[string("choose")][string("example")] = string(R"(n <- 5
k <- 2
x := choose(n, k))");
	help_strings[string("choose")][string("name")] = string(R"(choose)");
	help_arrays[string("clade")][string("authors")].push_back(string(R"(Will Pett)"));
	help_arrays[string("clade")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_arrays[string("clade")][string("authors")].push_back(string(R"(Michael Landis)"));
	help_strings[string("clade")][string("description")] = string(R"(A clade is a subtree within a phylogeny.)");
	help_strings[string("clade")][string("details")] = string(R"(Clades are defined in terms of a taxon set and a shared tree topology. In phylogenetic analyses, clades are generally used (a) to constrain tree topologies to match provided taxon relationships, (b) to identify the most recent common ancestor of a taxon set within a phylogeny, or (c) to apply node age calibrations on particular nodes in the phylogeny.)");
	help_strings[string("clade")][string("example")] = string(R"(# read in a tree
phy = readTrees("primates.tre")[1]
# get taxa from the tree
taxa = phy.taxa()
# create a clade for (1,2) using taxon objects
clade_12 = clade( taxa[1], taxa[2] )
# create a clade for (1,2,3)
clade_123 = clade( taxa[3], clade_12 )
# create a clade for (4,5) using taxon names
clade_45 = clade( "Taxon_4", "Taxon_5" )
# create a negative clade constraint
clade_not_23 = clade( taxa[2], taxa[3], negative=true )
# create an optional clade constraint
clade_67 = clade( taxa[6], taxa[7] )
clade_68 = clade( taxa[6], taxa[8] )
clade_67_or_68 = clade( clade_67, clade_68, optional_match=true ))");
	help_strings[string("clade")][string("name")] = string(R"(clade)");
	help_arrays[string("clade")][string("see_also")].push_back(string(R"(dnConstrainedTopology)"));
	help_arrays[string("clade")][string("see_also")].push_back(string(R"(tmrca)"));
	help_arrays[string("clade")][string("see_also")].push_back(string(R"(mrcaIndex)"));
	help_strings[string("clade")][string("title")] = string(R"(Clade)");
	help_arrays[string("clear")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("clear")][string("description")] = string(R"(Clear (e.g., remove) variables and functions from the workspace.)");
	help_strings[string("clear")][string("details")] = string(R"(The clear function removes either a given variable or all variables from the workspace. Clearing the workspace is very useful between analysis if you do not want to have old connections between variables hanging around.)");
	help_strings[string("clear")][string("example")] = string(R"(ls()   # check what is in the workspace
a <- 1
b := exp(a)
ls()   # check what is in the workspace
clear()
ls()   # check what is in the workspace
a <- 1
b := exp(a)
ls()   # check what is in the workspace
clear( b )
ls()   # check what is in the workspace)");
	help_strings[string("clear")][string("name")] = string(R"(clear)");
	help_arrays[string("clear")][string("see_also")].push_back(string(R"(exists)"));
	help_strings[string("clear")][string("title")] = string(R"(Clear the current workspace)");
	help_arrays[string("combineCharacter")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("combineCharacter")][string("description")] = string(R"(Creates a new data matrix by concatentating the provided data matrices (by order).)");
	help_strings[string("combineCharacter")][string("example")] = string(R"(# read in character data for locus_1
locus_1 = readContinuousCharacterData("locus_1.nex")
# read in character data for locus_2
locus_2 = readContinuousCharacterData("locus_2.nex")
# create concated locus for 1+2 (union of taxa)
locus_1_and_2 = concatenate( locus_1, locus_2 ))");
	help_strings[string("combineCharacter")][string("name")] = string(R"(combineCharacter)");
	help_strings[string("combineCharacter")][string("title")] = string(R"(Concatenate character matrices)");
	help_strings[string("computeWeightedNodeOrderConstraintsScore")][string("name")] = string(R"(computeWeightedNodeOrderConstraintsScore)");
	help_arrays[string("concatenate")][string("authors")].push_back(string(R"(Michael Landis)"));
	help_strings[string("concatenate")][string("description")] = string(R"(Creates a new data matrix by concatentating the provided data matrices (by order).)");
	help_strings[string("concatenate")][string("example")] = string(R"(# read in character data for locus_1
locus_1 = readDiscreteCharacterData("locus_1.nex")
# read in character data for locus_2
locus_2 = readDiscreteCharacterData("locus_2.nex")
# create concated locus for 1+2 (union of taxa)
locus_1_and_2 = concatenate( locus_1, locus_2 ))");
	help_strings[string("concatenate")][string("name")] = string(R"(concatenate)");
	help_strings[string("concatenate")][string("title")] = string(R"(Concatenate character matrices)");
	help_arrays[string("consensusTree")][string("authors")].push_back(string(R"(Seraina Klopfstein)"));
	help_arrays[string("consensusTree")][string("authors")].push_back(string(R"(Will Freyman)"));
	help_arrays[string("consensusTree")][string("authors")].push_back(string(R"(Will Pett)"));
	help_arrays[string("consensusTree")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("consensusTree")][string("description")] = string(R"(Calculates the majority-rule consensus topology from a trace of trees and summarizes branch lengths.)");
	help_strings[string("consensusTree")][string("example")] = string(R"(# Read in tree trace
tree_trace = readTreeTrace("output/my.trees", burnin=0.25)

# Generate the majority-rule consensus tree
map_tree = consensusTree(trace=tree_trace, cutoff=0.5, file="consensus.tree"))");
	help_strings[string("consensusTree")][string("name")] = string(R"(consensusTree)");
	help_arrays[string("consensusTree")][string("see_also")].push_back(string(R"(mapTree)"));
	help_arrays[string("consensusTree")][string("see_also")].push_back(string(R"(mccTree)"));
	help_arrays[string("consensusTree")][string("see_also")].push_back(string(R"(treeTrace)"));
	help_arrays[string("consensusTree")][string("see_also")].push_back(string(R"(readTreeTrace)"));
	help_strings[string("convertToPhylowood")][string("name")] = string(R"(convertToPhylowood)");
	help_strings[string("dfConstant")][string("name")] = string(R"(dfConstant)");
	help_strings[string("dfExponential")][string("name")] = string(R"(dfExponential)");
	help_strings[string("dfLinear")][string("name")] = string(R"(dfLinear)");
	help_arrays[string("diagonalMatrix")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("diagonalMatrix")][string("description")] = string(R"(Building a identity/diagonal matrix with 'n' columns and rows.)");
	help_strings[string("diagonalMatrix")][string("name")] = string(R"(diagonalMatrix)");
	help_arrays[string("dnBernoulli")][string("authors")].push_back(string(R"(John Huelsenbeck)"));
	help_strings[string("dnBernoulli")][string("description")] = string(R"(A Bernoulli-distributed random variable takes the value 1 with probability p and the value 0 with probability 1-p.)");
	help_strings[string("dnBernoulli")][string("example")] = string(R"(p ~ dnBeta(1.0,1.0)
x ~ dnBernoulli(p)
x.clamp(1)
moves[1] = mvSlide(p, delta=0.1, weight=1.0)
monitors[1] = screenmonitor(printgen=1000, separator = "        ", x)
mymodel = model(p)
mymcmc = mcmc(mymodel, monitors, moves)
mymcmc.burnin(generations=20000,tuningInterval=100)
mymcmc.run(generations=200000))");
	help_strings[string("dnBernoulli")][string("name")] = string(R"(dnBernoulli)");
	help_arrays[string("dnBernoulli")][string("see_also")].push_back(string(R"(dnBinomial)"));
	help_strings[string("dnBernoulli")][string("title")] = string(R"(Bernoulli Distribution)");
	help_arrays[string("dnBeta")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("dnBeta")][string("description")] = string(R"(The Beta probability distribution.)");
	help_strings[string("dnBeta")][string("details")] = string(R"(The Beta distribution is a probability distribution defined on 0 to 1.)");
	help_strings[string("dnBeta")][string("example")] = string(R"(p ~ dnBeta(1.0,1.0)
x ~ dnBernoulli(p)
x.clamp(1)
moves[1] = mvSlide(p, delta=0.1, weight=1.0)
monitors[1] = screenmonitor(printgen=1000, separator = "        ", x)
mymodel = model(p)
mymcmc = mcmc(mymodel, monitors, moves)
mymcmc.burnin(generations=20000,tuningInterval=100)
mymcmc.run(generations=200000))");
	help_strings[string("dnBeta")][string("name")] = string(R"(dnBeta)");
	help_arrays[string("dnBeta")][string("see_also")].push_back(string(R"(dnDirichlet)"));
	help_strings[string("dnBeta")][string("title")] = string(R"(Beta Distribution)");
	help_arrays[string("dnBimodalLognormal")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("dnBimodalLognormal")][string("description")] = string(R"(A bimodal lognormal distribution, that is, with probability p a value is distributed according to the first lognormal distribution and with probability 1-p from the second lognormal distribution.)");
	help_strings[string("dnBimodalLognormal")][string("example")] = string(R"(p ~ dnBeta(1.0,1.0)
x ~ dnBimodalLognormal(mean1=-1,mean2=1,sd1=0.1,sd2=0.1,p=p)
x.clamp( exp(1) )
moves[1] = mvSlide(p, delta=0.1, weight=1.0)
monitors[1] = screenmonitor(printgen=1000, separator = "        ", x)
mymodel = model(p)
mymcmc = mcmc(mymodel, monitors, moves)
mymcmc.burnin(generations=20000,tuningInterval=100)
mymcmc.run(generations=200000))");
	help_strings[string("dnBimodalLognormal")][string("name")] = string(R"(dnBimodalLognormal)");
	help_arrays[string("dnBimodalLognormal")][string("see_also")].push_back(string(R"(dnBimodalNormal)"));
	help_arrays[string("dnBimodalLognormal")][string("see_also")].push_back(string(R"(dnLognormal)"));
	help_strings[string("dnBimodalLognormal")][string("title")] = string(R"(Bimodal Lognormal Distribution)");
	help_arrays[string("dnBimodalNormal")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("dnBimodalNormal")][string("description")] = string(R"(A Bernoulli-distributed random variable takes the value 1 with probability p and the value 0 with probability 1-p.)");
	help_strings[string("dnBimodalNormal")][string("example")] = string(R"(p ~ dnBeta(1.0,1.0)
x ~ dnBimodalNormal(mean1=-1,mean2=1,sd1=0.1,sd2=0.1,p=p)
x.clamp( 1 )
moves[1] = mvSlide(p, delta=0.1, weight=1.0)
monitors[1] = screenmonitor(printgen=1000, separator = "        ", x)
mymodel = model(p)
mymcmc = mcmc(mymodel, monitors, moves)
mymcmc.burnin(generations=20000,tuningInterval=100)
mymcmc.run(generations=200000))");
	help_strings[string("dnBimodalNormal")][string("name")] = string(R"(dnBimodalNormal)");
	help_arrays[string("dnBimodalNormal")][string("see_also")].push_back(string(R"(dnBimodalLognormal)"));
	help_arrays[string("dnBimodalNormal")][string("see_also")].push_back(string(R"(dnNormal)"));
	help_strings[string("dnBimodalNormal")][string("title")] = string(R"(Bernoulli Distribution)");
	help_arrays[string("dnBinomial")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("dnBinomial")][string("description")] = string(R"(Binomial probability distribution of x successes in size trials.)");
	help_strings[string("dnBinomial")][string("details")] = string(R"(The binomial probability distribution defines the number of success in size trials, where each trial has the same success probability p. The probability is given by (size choose x) p^(x) * (1-p)^(size-p))");
	help_strings[string("dnBinomial")][string("example")] = string(R"(p ~ dnBeta(1.0,1.0)
x ~ dnBinomial(size=10,p)
x.clamp(8)
moves[1] = mvSlide(p, delta=0.1, weight=1.0)
monitors[1] = screenmonitor(printgen=1000, separator = "        ", x)
mymodel = model(p)
mymcmc = mcmc(mymodel, monitors, moves)
mymcmc.burnin(generations=20000,tuningInterval=100)
mymcmc.run(generations=200000))");
	help_strings[string("dnBinomial")][string("name")] = string(R"(dnBinomial)");
	help_arrays[string("dnBinomial")][string("see_also")].push_back(string(R"(dnBernoulli)"));
	help_strings[string("dnBinomial")][string("title")] = string(R"(Binomial Distribution)");
	help_strings[string("dnBirthDeath")][string("name")] = string(R"(dnBirthDeath)");
	help_strings[string("dnBirthDeathBurstProcess")][string("name")] = string(R"(dnBirthDeathBurstProcess)");
	help_strings[string("dnBirthDeathSamplingTreatment")][string("name")] = string(R"(dnBirthDeathSamplingTreatment)");
	help_arrays[string("dnBivariatePoisson")][string("authors")].push_back(string(R"(Alexander Zarebski)"));
	help_arrays[string("dnBivariatePoisson")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("dnBivariatePoisson")][string("description")] = string(R"(A Bivariate Poisson distribution defines probabilities for pairs of natural numbers.)");
	help_strings[string("dnBivariatePoisson")][string("example")] = string(R"(th1 ~ dnUniform(0.0,10.0)
th2 ~ dnUniform(0.0,10.0)
th0 ~ dnUniform(0.0,10.0)
 ~ dnBivariatePoisson(th1, th2, th0)
x.clamp([3, 3, 3])
oves[1] = mvSlide(th1, delta=0.01, weight=1.0)
moves[2] = mvSlide(th2, delta=0.01, weight=1.0)
oves[3] = mvSlide(th0, delta=0.01, weight=1.0)
monitors[1] = mnScreen(printgen=1000,  th0)
ymodel = model(th1)
mymcmc = mcmc(mymodel, monitors, moves)
ymcmc.burnin(generations=20000,tuningInterval=100)
mymcmc.run(generations=200000))");
	help_strings[string("dnBivariatePoisson")][string("name")] = string(R"(dnBivariatePoisson)");
	help_references[string("dnBivariatePoisson")].push_back(RbHelpReference(R"(Karlis D, Ntzoufras J (2003). Bayesian and Non-Bayesian Analysis of Soccer Data using Bivariate Poisson Regression Models.)",R"()",R"()"));
	help_arrays[string("dnBivariatePoisson")][string("see_also")].push_back(string(R"(dnPoisson)"));
	help_strings[string("dnBivariatePoisson")][string("title")] = string(R"(Bivariate Poisson Distribution)");
	help_strings[string("dnCBDSP")][string("name")] = string(R"(dnCBDSP)");
	help_strings[string("dnCDBDP")][string("name")] = string(R"(dnCDBDP)");
	help_arrays[string("dnCategorical")][string("authors")].push_back(string(R"(Fredrik Ronquist)"));
	help_strings[string("dnCategorical")][string("description")] = string(R"(The categorical distribution, sometimes referred to as the generalized Bernoulli distribution. It describes the probability of one of K different outcomes, labeled from 1 to K, with each outcome probability separately specified.)");
	help_strings[string("dnCategorical")][string("details")] = string(R"(The argument to the constructor is a simplex containing the probabilities of the outcomes. The outcomes are labeled from 1 to K, where K is the number of elements in the simplex. Outcome i has probability specified by component i in the simplex.

A typical scenario where a categorical variable is used is in the definition of a variable drawn from a mixture. A vector of mixture components is set up first, and then a stochastic variable drawn from a categorical distribution is used as an index in a deterministic assignment that points to a component in the mixture. See example below.)");
	help_strings[string("dnCategorical")][string("example")] = string(R"(# Define a stochastic variable x that is drawn from
# a categorical distribution with 4 categories, each
# category having the same probability, then examine
# the value of x.
x ~ dnCat( simplex(1,1,1,1) )
x

# Draw 10 values from the distribution and place them
# in a vector a, then examine a.
for ( i in 1:10 ) {
    a[i] <- x
    x.redraw()
}
a

# Use x in defining a deterministic variable y taking
# on values from a mixture of RealPos values representing
# rates from a discretized scaled gamma distribution
# with four categories.
shape ~ dnExp( 10.0 )
rates := fnDiscretizeGamma( shape, shape, 4 )
y := rates[x])");
	help_strings[string("dnCategorical")][string("name")] = string(R"(dnCategorical)");
	help_arrays[string("dnCategorical")][string("see_also")].push_back(string(R"(dnBinomial)"));
	help_strings[string("dnCategorical")][string("title")] = string(R"(The Categorical Distribution)");
	help_arrays[string("dnCauchy")][string("authors")].push_back(string(R"(Andrew Magee)"));
	help_strings[string("dnCauchy")][string("description")] = string(R"(Cauchy distribution with location equal to ‘location’ and scale equal to ‘scale’.)");
	help_strings[string("dnCauchy")][string("details")] = string(R"(The Cauchy distribution has density:

f(x) = 1/(pi * scale) * 1/(1 + x-(location/scale)^2))");
	help_strings[string("dnCauchy")][string("example")] = string(R"(# we simulate some obversations
x <- rCauchy(n=10,location=0,scale=1)
# let's see what the mean and the variance are.
The mean will not converge with more samples, the Cauchy family has no moments.
mean(x)
var(x)
sd(x))");
	help_strings[string("dnCauchy")][string("name")] = string(R"(dnCauchy)");
	help_strings[string("dnCauchy")][string("title")] = string(R"(Cauchy Distribution)");
	help_arrays[string("dnChisq")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("dnChisq")][string("description")] = string(R"(The chi-square probability distribution.)");
	help_strings[string("dnChisq")][string("example")] = string(R"(# The most important use of the chi-square distribution
# is arguable the quantile function.
# You can access it the following way:
df <- 10
a := qchisq(0.025, df)
a)");
	help_strings[string("dnChisq")][string("name")] = string(R"(dnChisq)");
	help_strings[string("dnChisq")][string("title")] = string(R"(Chi-Square Distribution)");
	help_strings[string("dnCoalescent")][string("name")] = string(R"(dnCoalescent)");
	help_strings[string("dnCoalescentDemography")][string("name")] = string(R"(dnCoalescentDemography)");
	help_strings[string("dnCoalescentSkyline")][string("name")] = string(R"(dnCoalescentSkyline)");
	help_strings[string("dnCompleteBirthDeath")][string("name")] = string(R"(dnCompleteBirthDeath)");
	help_strings[string("dnConstrainedNodeAge")][string("name")] = string(R"(dnConstrainedNodeAge)");
	help_strings[string("dnConstrainedNodeOrder")][string("name")] = string(R"(dnConstrainedNodeOrder)");
	help_strings[string("dnConstrainedTopology")][string("name")] = string(R"(dnConstrainedTopology)");
	help_strings[string("dnCppNormal")][string("name")] = string(R"(dnCppNormal)");
	help_strings[string("dnDPP")][string("name")] = string(R"(dnDPP)");
	help_strings[string("dnDecomposedInvWishart")][string("name")] = string(R"(dnDecomposedInvWishart)");
	help_arrays[string("dnDirichlet")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("dnDirichlet")][string("description")] = string(R"(Dirichlet probability distribution on a simplex.)");
	help_strings[string("dnDirichlet")][string("details")] = string(R"(The Dirichlet probability distribution is the generalization of the Beta distribution. A random variable from a Dirichlet distribution is a simplex, i.e., a vector of probabilities that sum to 1.)");
	help_strings[string("dnDirichlet")][string("example")] = string(R"(# lets get a draw from a Dirichlet distribution
a <- [1,1,1,1]   # we could also use rep(1,4)
b ~ dnDirichlet(a)
b
# let check if b really sums to 1
sum(b))");
	help_strings[string("dnDirichlet")][string("name")] = string(R"(dnDirichlet)");
	help_arrays[string("dnDirichlet")][string("see_also")].push_back(string(R"(simplex)"));
	help_strings[string("dnDirichlet")][string("title")] = string(R"(Dirichlet Distribution)");
	help_strings[string("dnDiversityDependentYule")][string("name")] = string(R"(dnDiversityDependentYule)");
	help_arrays[string("dnDuplicationLoss")][string("authors")].push_back(string(R"(Sebastian Hoehna, Bastien Boussau, Dominik ...)"));
	help_strings[string("dnDuplicationLoss")][string("description")] = string(R"(Multispecies coalescent distribution describing how gene trees can be generated from within a species tree given a constant effective population size. Requires an ultrametric species tree, a single effective population size (a single real positive), and taxa with species and individual names.)");
	help_strings[string("dnDuplicationLoss")][string("details")] = string(R"(The species tree must be ultrametric.
The effective population size is constant across the species tree.)");
	help_strings[string("dnDuplicationLoss")][string("example")] = string(R"(# We are going to save the trees we simulate in the folder simulatedTrees:
dataFolder = "simulatedTrees/"
# Let’s simulate a species tree with 10 taxa, 2 gene trees, 3 alleles per species:
n_species <- 10
n_genes <- 2
n_alleles <- 3
# we simulate an ultrametric species tree:
# Species names:
for (i in 1:n_species) {
        species[i] <- taxon(taxonName="Species_"+i, speciesName="Species_"+i)
}
spTree ~ dnBirthDeath(lambda=0.3, mu=0.2, rootAge=10, rho=1, samplingStrategy="uniform", condition="nTaxa", taxa=species)
print(spTree)
# let's pick a constant effective population size of 50:
popSize <- 50
# let's simulate gene trees now:
# taxa names:
for (g in 1:n_genes) {
  for (i in 1:n_species) {
    for (j in 1:n_alleles) {
        taxons[g][(i-1)*n_alleles+j] <- taxon(taxonName="Species_"+i+"_"+j, speciesName="Species_"+i)
    }
  }
  geneTrees[g] ~ dnMultiSpeciesCoalescent(speciesTree=spTree, Ne=popSize, taxa=taxons[g])
  print(geneTrees[g])
}
# We can save the species tree and the gene trees:
write(spTree, filename=dataFolder+"speciesTree")
# Saving the gene trees
for (i in 1:(n_genes)) {
  write(geneTrees[i], filename=dataFolder+"geneTree_"+i+".tree")
})");
	help_strings[string("dnDuplicationLoss")][string("name")] = string(R"(dnDuplicationLoss)");
	help_references[string("dnDuplicationLoss")].push_back(RbHelpReference(R"(Bayes Estimation of Species Divergence Times and Ancestral Population Sizes Using DNA Sequences From Multiple Loci. Bruce Rannala and Ziheng Yang. GENETICS August 1, 2003 vol. 164 no. 4 1645-1656.)",R"()",R"(http://www.genetics.org/content/164/4/1645.short )"));
	help_arrays[string("dnDuplicationLoss")][string("see_also")].push_back(string(R"(dnMultiSpeciesCoalescentUniformPrior)"));
	help_arrays[string("dnDuplicationLoss")][string("see_also")].push_back(string(R"(dnMultiSpeciesCoalescentInverseGamma)"));
	help_strings[string("dnDuplicationLoss")][string("title")] = string(R"(Multispecies coalescent Distribution)");
	help_strings[string("dnEmpiricalSample")][string("name")] = string(R"(dnEmpiricalSample)");
	help_arrays[string("dnEmpiricalTree")][string("authors")].push_back(string(R"(Will Freyman)"));
	help_arrays[string("dnEmpiricalTree")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_arrays[string("dnEmpiricalTree")][string("authors")].push_back(string(R"(Will Pett)"));
	help_strings[string("dnEmpiricalTree")][string("description")] = string(R"(Creates a distribution of trees from a trace of trees.)");
	help_strings[string("dnEmpiricalTree")][string("example")] = string(R"(# Read in tree trace
tree_trace = readTreeTrace("output/my.trees", burnin=0.25)

# Create a distribution of trees
tree ~ dnEmpiricalTree(tree_trace)

# Add an MCMC move
moves[1] = mvEmpiricalTree(tree))");
	help_strings[string("dnEmpiricalTree")][string("name")] = string(R"(dnEmpiricalTree)");
	help_arrays[string("dnEmpiricalTree")][string("see_also")].push_back(string(R"(mvEmpiricalTree)"));
	help_arrays[string("dnEmpiricalTree")][string("see_also")].push_back(string(R"(treeTrace)"));
	help_arrays[string("dnEmpiricalTree")][string("see_also")].push_back(string(R"(readTreeTrace)"));
	help_strings[string("dnEpisodicBirthDeath")][string("name")] = string(R"(dnEpisodicBirthDeath)");
	help_strings[string("dnEvent")][string("name")] = string(R"(dnEvent)");
	help_arrays[string("dnExponential")][string("authors")].push_back(string(R"(Michael Landis)"));
	help_strings[string("dnExponential")][string("description")] = string(R"(Exponential distribution with rate equal to ‘lambda’.)");
	help_strings[string("dnExponential")][string("details")] = string(R"(The exponential distribution has density:

f(x) = lambda * exp(-lambda*x)

where lambda is the rate parameter.)");
	help_strings[string("dnExponential")][string("example")] = string(R"(# we set a rate parameter
rate <- 10.0
# we create an exponentially distributed random variable
x ~ dnExponential(lambda=rate)
# compute the probability of the variable
x.probability())");
	help_strings[string("dnExponential")][string("name")] = string(R"(dnExponential)");
	help_strings[string("dnExponential")][string("title")] = string(R"(Exponential Distribution)");
	help_arrays[string("dnFossilizedBirthDeathRange")][string("authors")].push_back(string(R"(Walker Pett)"));
	help_strings[string("dnFossilizedBirthDeathRange")][string("description")] = string(R"(The fossilized birth death range process (FBDRP) describes the distribution of a matrix of species origination and extinction times under a model of asymmetric speciation and sampling of extinct species.)");
	help_strings[string("dnFossilizedBirthDeathRange")][string("details")] = string(R"(Fossil species are represented by a collection of fossil occurrences with uncertainty. Speciation, extinction and sampling rates may be time-homogeneous or piecewise time-heterogeneous. If time-heterogeneous rates are provided, then a vector of rate change time-points musts also be provided. If only a subset of fossil occurrences is provided for each species (e.g. only first and last occurrencces), then the remaining number of fossil samples may be marginalized by specifying `complete=FALSE`. Under the hood, the fossil data is augmented with oldest occurrence ages for each species, which must be sampled during MCMC using `mvResampleFBD`. Setting `BDS` to true causes the model to assume complete lineage sampling i.e. using the Birth-Death with Rateshifts (BDS) model of Silvestro et al. (2019).)");
	help_strings[string("dnFossilizedBirthDeathRange")][string("example")] = string(R"(lambda ~ dnExp(10)
mu ~ dnExp(10)
psi ~ dnExp(10)

bd ~ dnFBDRP(lambda=lambda, mu=mu, psi=psi, rho=1, taxa=taxa)

moves.append( mvMatrixElementScale(bd, weight=taxa.size()) )
moves.append( mvMatrixElementSlide(bd, weight=taxa.size()) ))");
	help_strings[string("dnFossilizedBirthDeathRange")][string("name")] = string(R"(dnFossilizedBirthDeathRange)");
	help_references[string("dnFossilizedBirthDeathRange")].push_back(RbHelpReference(R"(The fossilized birth-death model for the analysis of stratigraphic range data under different speciation modes. Stadler, Tanja et al. Journal of theoretical biology, 447:41-55.)",R"()",R"(https://www.sciencedirect.com/science/article/pii/S002251931830119X )"));
	help_references[string("dnFossilizedBirthDeathRange")].push_back(RbHelpReference(R"(Improved estimation of macroevolutionary rates from fossil data using a Bayesian framework. Silvestro, Daniele et al. Paleobiology, 45:546-570.)",R"(https://doi.org/10.1017/pab.2019.23)",R"(https://www.cambridge.org/core/journals/paleobiology/article/improved-estimation-of-macroevolutionary-rates-from-fossil-data-using-a-bayesian-framework/334F08A74A6C92F1FEAD91A71FE59A1C )"));
	help_arrays[string("dnFossilizedBirthDeathRange")][string("see_also")].push_back(string(R"(dnBirthDeathSamplingTreatment)"));
	help_arrays[string("dnFossilizedBirthDeathRange")][string("see_also")].push_back(string(R"(mvResampleFBD)"));
	help_arrays[string("dnFossilizedBirthDeathSpeciation")][string("authors")].push_back(string(R"(Walker Pett)"));
	help_strings[string("dnFossilizedBirthDeathSpeciation")][string("description")] = string(R"(The fossilized birth death speciation process (FBDSP) describes the diversification and sampling of extant and extinct species trees under a mixed model of asymmetric, symmetric and anagenetic speciation.)");
	help_strings[string("dnFossilizedBirthDeathSpeciation")][string("details")] = string(R"(Fossil species are represented by a collection of fossil occurrences with uncertainty. Speciation, extinction and sampling rates may be time-homogeneous or piecewise time-heterogeneous. If time-heterogeneous rates are provided, then a vector of rate change time-points musts also be provided. If only a subset of fossil occurrences is provided for each species (e.g. only first and last occurrencces), then the remaining number of fossil samples may be marginalized by specifying `complete=FALSE`. Under the hood, the fossil data is augmented with oldest occurrence ages for each species, which must be sampled during MCMC using `mvResampleFBD`. Tips represent extinction events, and therefore should be sampled during MCMC using e.g. `mvTipTimeSlideUniform`.)");
	help_strings[string("dnFossilizedBirthDeathSpeciation")][string("example")] = string(R"(lambda ~ dnExp(10)
mu ~ dnExp(10)
psi ~ dnExp(10)

min_age = 0.0
for(i in 1:taxa.size())
{
        if ( taxa[i].getMinAge() > min_age )
        {
                min_age = taxa[i].getMinAge()
        }
}

origin ~ dnExp(1/10)

bd ~ dnFBDSP(originAge=min_age+origin, lambda=lambda, mu=mu, psi=psi, rho=1, taxa=taxa, complete=FALSE)

moves.append( mvFNPR(bd, weight = taxa.size()) )
moves.append( mvNodeTimeSlideUniform(bd, weight = taxa.size()) )
moves.append( mvRootTimeSlideUniform(bd, origin=origin, weight = taxa.size()) )
moves.append( mvTipTimeSlideUniform(bd, weight = taxa.size()) ))");
	help_strings[string("dnFossilizedBirthDeathSpeciation")][string("name")] = string(R"(dnFossilizedBirthDeathSpeciation)");
	help_references[string("dnFossilizedBirthDeathSpeciation")].push_back(RbHelpReference(R"(The fossilized birth-death model for the analysis of stratigraphic range data under different speciation modes. Stadler, Tanja et al. Journal of theoretical biology, 447:41-55.)",R"()",R"(https://www.sciencedirect.com/science/article/pii/S002251931830119X )"));
	help_arrays[string("dnFossilizedBirthDeathSpeciation")][string("see_also")].push_back(string(R"(dnFossilizedBirthDeathRange)"));
	help_arrays[string("dnFossilizedBirthDeathSpeciation")][string("see_also")].push_back(string(R"(dnBirthDeathSamplingTreatment)"));
	help_arrays[string("dnFossilizedBirthDeathSpeciation")][string("see_also")].push_back(string(R"(mvResampleFBD)"));
	help_arrays[string("dnGamma")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("dnGamma")][string("description")] = string(R"(Gamma probability distribution for positive real numbers.)");
	help_strings[string("dnGamma")][string("details")] = string(R"(The gamma distribution is the probability of the sum of exponentially distributed variables. Thus, it provides a natural prior distribution for parameters that could be considered as sums of exponential variables.)");
	help_strings[string("dnGamma")][string("example")] = string(R"(# lets simulate
a <- rgamma(1000,shape=4,rate=4)
# we expect a mean of 1
mean(a)

# create a random variable
x ~ dnGamma(shape=4,rate=1)
x)");
	help_strings[string("dnGamma")][string("name")] = string(R"(dnGamma)");
	help_arrays[string("dnGamma")][string("see_also")].push_back(string(R"(dnExponential)"));
	help_strings[string("dnGamma")][string("title")] = string(R"(Gamma Distribution)");
	help_arrays[string("dnGeometric")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("dnGeometric")][string("description")] = string(R"(A geometric distribution defines a random variable on natural numbers. The geometric distribution describes the number of success until the first failure, with success probability p.)");
	help_strings[string("dnGeometric")][string("example")] = string(R"(p ~ dnBeta(1.0,1.0)
x ~ dnGeom(p)
x.clamp(10)
moves[1] = mvSlide(p, delta=0.1, weight=1.0)
monitors[1] = screenmonitor(printgen=1000, separator = "        ", p)
mymodel = model(p)
mymcmc = mcmc(mymodel, monitors, moves)
mymcmc.burnin(generations=20000,tuningInterval=100)
mymcmc.run(generations=200000))");
	help_strings[string("dnGeometric")][string("name")] = string(R"(dnGeometric)");
	help_arrays[string("dnGeometric")][string("see_also")].push_back(string(R"(dnPoisson)"));
	help_arrays[string("dnGeometric")][string("see_also")].push_back(string(R"(mvRandomIntegerWalk)"));
	help_strings[string("dnGeometric")][string("title")] = string(R"(Geometric Distribution)");
	help_strings[string("dnGilbertGraph")][string("name")] = string(R"(dnGilbertGraph)");
	help_strings[string("dnHBDP")][string("name")] = string(R"(dnHBDP)");
	help_arrays[string("dnHalfCauchy")][string("authors")].push_back(string(R"(Andrew Magee)"));
	help_strings[string("dnHalfCauchy")][string("description")] = string(R"(Half-Cauchy distribution with location equal to ‘location’ and scale equal to ‘scale’.)");
	help_strings[string("dnHalfCauchy")][string("details")] = string(R"(The half-Cauchy distribution has density:

f(x) = 2/(pi * sigma) * 1/(1 + x-(location/scale)^2))");
	help_strings[string("dnHalfCauchy")][string("example")] = string(R"(# we simulate some obversations
x <- rHalfCauchy(n=10,location=0,scale=1)
# let's see what the minimum is (you could do the max too). If this is not ‘location’, we're in trouble!
min(x)
# let's also see what the mean and the variance are.
The mean will not converge with more samples, the Cauchy family has no moments.
mean(x)
var(x)
sd(x))");
	help_strings[string("dnHalfCauchy")][string("name")] = string(R"(dnHalfCauchy)");
	help_strings[string("dnHalfCauchy")][string("title")] = string(R"(half-Cauchy Distribution)");
	help_arrays[string("dnHalfNormal")][string("authors")].push_back(string(R"(Andrew Magee)"));
	help_strings[string("dnHalfNormal")][string("description")] = string(R"(half-normal (gaussian) distribution with offset equal to ‘offset’ and standard deviation equal to ‘sd’.)");
	help_strings[string("dnHalfNormal")][string("details")] = string(R"(The half-normal distribution has density:

 f(x) = 2/(sqrt(2 pi) sigma) e^-((x - offset)^2/(2 sigma^2)) where offset is the offset of the distribution and sigma the standard deviation.

f(x) = 2/(sqrt(2 pi) sigma) e^-((x - offset)^2/(2 sigma^2))

where offset is the offset of the distribution and sigma the standard deviation.)");
	help_strings[string("dnHalfNormal")][string("example")] = string(R"(# we simulate some oversations
x <- rhalfNormal(n=10,offset=0,sd=10)
# let's see what the minimum is (you could do the max too)
# the minimum should never be less than the offset
min(x)
# let's also see what the mean and the variance are
mean(x)
var(x)
sd(x))");
	help_strings[string("dnHalfNormal")][string("name")] = string(R"(dnHalfNormal)");
	help_arrays[string("dnHalfNormal")][string("see_also")].push_back(string(R"(dnNormal)"));
	help_arrays[string("dnHalfNormal")][string("see_also")].push_back(string(R"(dnLognormal)"));
	help_strings[string("dnHalfNormal")][string("title")] = string(R"(half-Normal Distribution)");
	help_strings[string("dnHeterochronousCoalescent")][string("name")] = string(R"(dnHeterochronousCoalescent)");
	help_strings[string("dnHeterochronousCoalescentSkyline")][string("name")] = string(R"(dnHeterochronousCoalescentSkyline)");
	help_arrays[string("dnInverseGamma")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("dnInverseGamma")][string("description")] = string(R"(inverse-gamma probability distribution for positive real numbers.)");
	help_strings[string("dnInverseGamma")][string("details")] = string(R"(The inverse Gamma distribution is the probability of the sum of exponentially distributed variables. Thus, it provides a natural prior distribution for parameters that could be considered as sums of exponential variables.)");
	help_strings[string("dnInverseGamma")][string("example")] = string(R"(# lets simulate
a <- rinverseGamma(1000,shape=4,rate=4)
# we expect a mean of 1
mean(a)

# create a random variable
x ~ dnInverseGamma(shape=4,rate=1)
x)");
	help_strings[string("dnInverseGamma")][string("name")] = string(R"(dnInverseGamma)");
	help_arrays[string("dnInverseGamma")][string("see_also")].push_back(string(R"(dnExponential)"));
	help_strings[string("dnInverseGamma")][string("title")] = string(R"(inverseGamma Distribution)");
	help_strings[string("dnInverseWishart")][string("name")] = string(R"(dnInverseWishart)");
	help_arrays[string("dnLKJ")][string("authors")].push_back(string(R"(Michael R. May)"));
	help_strings[string("dnLKJ")][string("description")] = string(R"(The LKJ (Lewandowski-Kurowicka-Joe) distribution on correlation matrices with concentration parameter eta.)");
	help_strings[string("dnLKJ")][string("details")] = string(R"(The LKJ distribution is uniform over positive-definite correlation matrices when eta=1.The probability density of a correlation matrix under the LKJ distribution is:f(x) = det(x)^(eta - 1))");
	help_strings[string("dnLKJ")][string("example")] = string(R"(
# we simulate a correlation matrix.
R <- rLKJ(n=1, eta=1, dim=5)

# let's print the simulated correlation matrix
R)");
	help_strings[string("dnLKJ")][string("name")] = string(R"(dnLKJ)");
	help_references[string("dnLKJ")].push_back(RbHelpReference(R"(Lewandowski D, D Kurowicka, H Joe (2009). Generating random correlation matrices based on vines and extended onion method. Journal of Multivariate Analysis, 100(9):1989-2001.)",R"()",R"()"));
	help_arrays[string("dnLKJ")][string("see_also")].push_back(string(R"(dnLKJPartial)"));
	help_strings[string("dnLKJ")][string("title")] = string(R"(LKJ Distribution)");
	help_arrays[string("dnLKJPartial")][string("authors")].push_back(string(R"(Michael R. May)"));
	help_strings[string("dnLKJPartial")][string("description")] = string(R"(The LKJ (Lewandowski-Kurowicka-Joe) distribution (on the partial correlation matrix) with concentration parameter eta.)");
	help_strings[string("dnLKJPartial")][string("details")] = string(R"(The LKJPartial distribution is uniform over positive-definite correlation matrices when eta=1.

The probability density of a correlation matrix under the LKJ distribution is:

f(x) = det(x)^(eta - 1))");
	help_strings[string("dnLKJPartial")][string("example")] = string(R"(# we simulate a partial correlation matrix.
P <- rLKJPartial(n=1, eta=1, dim=5)

# let's print the simulated partial correlation matrix
P)");
	help_strings[string("dnLKJPartial")][string("name")] = string(R"(dnLKJPartial)");
	help_references[string("dnLKJPartial")].push_back(RbHelpReference(R"(Lewandowski D, D Kurowicka, H Joe (2009). Generating random correlation matrices based on vines and extended onion method. Journal of Multivariate Analysis, 100(9):1989-2001.)",R"()",R"()"));
	help_arrays[string("dnLKJPartial")][string("see_also")].push_back(string(R"(dnLKJ)"));
	help_strings[string("dnLKJPartial")][string("title")] = string(R"(LKJ Distribution (for partial correlation matrices))");
	help_arrays[string("dnLaplace")][string("authors")].push_back(string(R"(Will Freyman)"));
	help_strings[string("dnLaplace")][string("description")] = string(R"(Laplace distribution with mean equal to ‘mean’ and scale equal to ‘scale’.)");
	help_strings[string("dnLaplace")][string("details")] = string(R"(The Laplace distribution has density:

f(x) = 1/(2 b) e^-(abs(x-mu)/b)

where mu is the mean of the distribution and b the scale.)");
	help_strings[string("dnLaplace")][string("name")] = string(R"(dnLaplace)");
	help_arrays[string("dnLaplace")][string("see_also")].push_back(string(R"(dnExponential)"));
	help_arrays[string("dnLaplace")][string("see_also")].push_back(string(R"(dnNormal)"));
	help_strings[string("dnLaplace")][string("title")] = string(R"(Laplace Distribution)");
	help_arrays[string("dnLogExponential")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("dnLogExponential")][string("description")] = string(R"(A real number x has a log-Exponential distribution if y = exp(x) has Exponential distribution.)");
	help_strings[string("dnLogExponential")][string("details")] = string(R"(The log-Exponential distribution is defined over real numbers. Saying that x is log-Exponential is equivalent to saying that y = exp(x) is Exponential. The log-Exponential distribution therefore expresses lack of information about the order of magnitude of a scale parameter:  if x has a log-Exponential distribution, then it has equal chance to be contained by any of the intervals of the form (10^k, 10^(k+1)) within the allowed range.

The density is p(x) = ???, which can be seen by defining x = ln(y) where y has Exponential distribution and apply the change-of-variable formula.)");
	help_strings[string("dnLogExponential")][string("example")] = string(R"(# a log-Exponential prior over the rate of change of a Brownian trait (or a Brownian relaxed clock)
trueTree = readTrees("data/primates.tree")[1]
log_sigma ~ dnLogExponential(lambda=1)
sigma := exp(log_sigma)
X ~ dnBrownian(trueTree,sigma)
# ...)");
	help_strings[string("dnLogExponential")][string("name")] = string(R"(dnLogExponential)");
	help_arrays[string("dnLogExponential")][string("see_also")].push_back(string(R"(dnExponential)"));
	help_strings[string("dnLogExponential")][string("title")] = string(R"(Log-Exponential Distribution)");
	help_arrays[string("dnLognormal")][string("authors")].push_back(string(R"(Michael Landis)"));
	help_strings[string("dnLognormal")][string("description")] = string(R"(Lognormal distribution is the distribution for a log-transformed normally distributed random variable with mean 'mu' and standard deviation 'sigma'.)");
	help_strings[string("dnLognormal")][string("details")] = string(R"(The lognormal random variable is defined as

:X = exp(mu + sigma Z)

where mu is the mean parameter, sigma is the standard deviation, and Z is a standard normal random variable. Note, in effect, the mean and standard deviation provide the location and scale of the exponentiated normal variate, mu + sigma Z.The lognormal distribution has density:

f(x) = 1/(x sigma sqrt(2 pi)) e^-((ln x - mu)^2/(2 sigma^2))

where mu is the mean of the distribution and sigma the standard deviation.)");
	help_strings[string("dnLognormal")][string("example")] = string(R"(# set an expected value for x
expectation_of_x <- 1
# set a mean and sd parameter
sd <- 0.5
mean <- ln(expectation_of_x) - 0.5 * sd^2
# create a lognormal distribution with expected value of 1
x ~ dnLognormal(mean=mean, sd=sd))");
	help_strings[string("dnLognormal")][string("name")] = string(R"(dnLognormal)");
	help_strings[string("dnLognormal")][string("title")] = string(R"(Lognormal Distribution)");
	help_arrays[string("dnLoguniform")][string("authors")].push_back(string(R"(Nicolas Lartillot)"));
	help_strings[string("dnLoguniform")][string("description")] = string(R"(A strictly positive real number x has a log-uniform distribution over interval (min,max) if its logarithm y = ln(x) has uniform distribution over interval (ln(min),ln(max)).)");
	help_strings[string("dnLoguniform")][string("details")] = string(R"(The log-uniform distribution is defined over strictly positive real numbers. Saying that x is log-uniform is equivalent to saying that y = ln(x) is uniform. The log-uniform distribution therefore expresses lack of information about the order of magnitude of a scale parameter:  if x has a log-uniform distribution, then it has equal chance to be contained by any of the intervals of the form (10^k, 10^(k+1)) within the allowed range.

The density is p(x) = 1/x, which can be seen by defining x = exp(y) where y has uniform distribution and apply the change-of-variable formula.

The log-uniform distribution is improper when defined over the entire positive real line. To always make it proper, in RevBayes, a min and a max should always be specified.)");
	help_strings[string("dnLoguniform")][string("example")] = string(R"(# a log-uniform prior over the rate of change of a Brownian trait (or a Brownian relaxed clock)
trueTree = readTrees("data/primates.tree")[1]
sigma ~ dnLogUniform(min=0.001, max=1000)
X ~ dnBrownian(trueTree,sigma)
# ...)");
	help_strings[string("dnLoguniform")][string("name")] = string(R"(dnLoguniform)");
	help_arrays[string("dnLoguniform")][string("see_also")].push_back(string(R"(dnUniform)"));
	help_strings[string("dnLoguniform")][string("title")] = string(R"(Log-Uniform Distribution)");
	help_strings[string("dnMixture")][string("name")] = string(R"(dnMixture)");
	help_strings[string("dnMixtureVector")][string("name")] = string(R"(dnMixtureVector)");
	help_arrays[string("dnMultiSpeciesCoalescent")][string("authors")].push_back(string(R"(Sebastian Hoehna, Bastien Boussau)"));
	help_strings[string("dnMultiSpeciesCoalescent")][string("description")] = string(R"(Multispecies coalescent distribution describing how gene trees can be generated from within a species tree given a constant effective population size. Requires an ultrametric species tree, a single effective population size (a single real positive), and taxa with species and individual names.)");
	help_strings[string("dnMultiSpeciesCoalescent")][string("details")] = string(R"(The species tree must be ultrametric.
The effective population size is constant across the species tree.)");
	help_strings[string("dnMultiSpeciesCoalescent")][string("example")] = string(R"(# We are going to save the trees we simulate in the folder simulatedTrees:
dataFolder = "simulatedTrees/"
# Let’s simulate a species tree with 10 taxa, 2 gene trees, 3 alleles per species:
n_species <- 10
n_genes <- 2
n_alleles <- 3
# we simulate an ultrametric species tree:
# Species names:
for (i in 1:n_species) {
        species[i] <- taxon(taxonName="Species_"+i, speciesName="Species_"+i)
}
spTree ~ dnBirthDeath(lambda=0.3, mu=0.2, rootAge=10, rho=1, samplingStrategy="uniform", condition="nTaxa", taxa=species)
print(spTree)
# let's pick a constant effective population size of 50:
popSize <- 50
# let's simulate gene trees now:
# taxa names:
for (g in 1:n_genes) {
  for (i in 1:n_species) {
    for (j in 1:n_alleles) {
        taxons[g][(i-1)*n_alleles+j] <- taxon(taxonName="Species_"+i+"_"+j, speciesName="Species_"+i)
    }
  }
  geneTrees[g] ~ dnMultiSpeciesCoalescent(speciesTree=spTree, Ne=popSize, taxa=taxons[g])
  print(geneTrees[g])
}
# We can save the species tree and the gene trees:
write(spTree, filename=dataFolder+"speciesTree")
# Saving the gene trees
for (i in 1:(n_genes)) {
  write(geneTrees[i], filename=dataFolder+"geneTree_"+i+".tree")
})");
	help_strings[string("dnMultiSpeciesCoalescent")][string("name")] = string(R"(dnMultiSpeciesCoalescent)");
	help_references[string("dnMultiSpeciesCoalescent")].push_back(RbHelpReference(R"(Bayes Estimation of Species Divergence Times and Ancestral Population Sizes Using DNA Sequences From Multiple Loci. Bruce Rannala and Ziheng Yang. GENETICS August 1, 2003 vol. 164 no. 4 1645-1656.)",R"()",R"(http://www.genetics.org/content/164/4/1645.short )"));
	help_arrays[string("dnMultiSpeciesCoalescent")][string("see_also")].push_back(string(R"(dnMultiSpeciesCoalescentUniformPrior)"));
	help_arrays[string("dnMultiSpeciesCoalescent")][string("see_also")].push_back(string(R"(dnMultiSpeciesCoalescentInverseGamma)"));
	help_strings[string("dnMultiSpeciesCoalescent")][string("title")] = string(R"(Multispecies coalescent Distribution)");
	help_arrays[string("dnMultiSpeciesCoalescentInverseGamma")][string("authors")].push_back(string(R"(Sebastian Hoehna, Bastien Boussau)"));
	help_strings[string("dnMultiSpeciesCoalescentInverseGamma")][string("description")] = string(R"(Multispecies coalescent distribution describing how gene trees can be generated from within a species tree given effective population sizes. Requires an ultrametric species tree, parameters of an inverse gamma prior on effective population sizes, and taxa with species and individual names.)");
	help_strings[string("dnMultiSpeciesCoalescentInverseGamma")][string("details")] = string(R"(The species tree must be ultrametric.
Parameters of an inverse gamma prior on effective population sizes must be provided.
This distribution uses a conjugate prior on effective population sizes. As a consequence, effective population sizes are integrated out and treated as nuisance parameters.

If you are interested in reconstructing ancestral effective population sizes, use dnMultiSpeciesCoalescent.)");
	help_strings[string("dnMultiSpeciesCoalescentInverseGamma")][string("example")] = string(R"(# We are going to save the trees we simulate in the folder simulatedTrees:
dataFolder = "simulatedTrees/"
# Let’s simulate a species tree with 10 taxa, 2 gene trees, 3 alleles per species:
n_species <- 10
n_genes <- 2
n_alleles <- 3
# we simulate an ultrametric species tree:
# Species names:
for (i in 1:n_species) {
        species[i] <- taxon(taxonName="Species_"+i, speciesName="Species_"+i)
}
spTree ~ dnBirthDeath(lambda=0.3, mu=0.2, rootAge=10, rho=1, samplingStrategy="uniform", condition="nTaxa", taxa=species)
print(spTree)
# let's pick constant parameters for the inverse gamma distribution:
alpha <- 3
beta <- 0.003
# let's simulate gene trees now:
# taxa names:
for (g in 1:n_genes) {
  for (i in 1:n_species) {
    for (j in 1:n_alleles) {
        taxons[g][(i-1)*n_alleles+j] <- taxon(taxonName="Species_"+i+"_"+j, speciesName="Species_"+i)
    }
  }
  geneTrees[g] ~ dnMultiSpeciesCoalescentInverseGamma(speciesTree=spTree, shape=alpha, scale=beta, taxa=taxons[g])
  print(geneTrees[g])
}
# We can save the species tree and the gene trees:
write(spTree, filename=dataFolder+"speciesTree")
# Saving the gene trees
for (i in 1:(n_genes)) {
  write(geneTrees[i], filename=dataFolder+"geneTree_"+i+".tree")
})");
	help_strings[string("dnMultiSpeciesCoalescentInverseGamma")][string("name")] = string(R"(dnMultiSpeciesCoalescentInverseGamma)");
	help_references[string("dnMultiSpeciesCoalescentInverseGamma")].push_back(RbHelpReference(R"(' Algorithmic improvements to species delimitation and phylogeny estimation under the multispecies coalescent. Jones G.  Journal of Mathematical Biology. 2016.')",R"('DOI: 10.1007/s00285-016-1034-0')",R"(http://www.indriid.com/2016/2016-06-01-STACEY.pdf )"));
	help_arrays[string("dnMultiSpeciesCoalescentInverseGamma")][string("see_also")].push_back(string(R"(dnMultiSpeciesCoalescent)"));
	help_arrays[string("dnMultiSpeciesCoalescentInverseGamma")][string("see_also")].push_back(string(R"(dnMultiSpeciesCoalescentUniformPrior)"));
	help_strings[string("dnMultiSpeciesCoalescentInverseGamma")][string("title")] = string(R"(Multispecies coalescent Distribution with inverse gamma prior on effective population sizes)");
	help_arrays[string("dnMultiSpeciesCoalescentUniformPrior")][string("authors")].push_back(string(R"(Sebastian Hoehna, Bastien Boussau)"));
	help_strings[string("dnMultiSpeciesCoalescentUniformPrior")][string("description")] = string(R"(Multispecies coalescent distribution describing how gene trees can be generated from within a species tree given effective population sizes. Requires an ultrametric species tree, effective population size(s) (a single real positive or a vector of real positives), and taxa with species and individual names.)");
	help_strings[string("dnMultiSpeciesCoalescentUniformPrior")][string("details")] = string(R"(The species tree must be ultrametric.
Effective population sizes can be constant across the species tree, if a single real positive is provided, or branchwise, if a vector is provided.)");
	help_strings[string("dnMultiSpeciesCoalescentUniformPrior")][string("example")] = string(R"(# We are going to save the trees we simulate in the folder simulatedTrees:
dataFolder = "simulatedTrees/"
# Let’s simulate a species tree with 10 taxa, 2 gene trees, 3 alleles per species:
n_species <- 10
n_genes <- 2
n_alleles <- 3
# we simulate an ultrametric species tree:
# Species names:
for (i in 1:n_species) {
        species[i] <- taxon(taxonName="Species_"+i, speciesName="Species_"+i)
}
spTree ~ dnBirthDeath(lambda=0.3, mu=0.2, rootAge=10, rho=1, samplingStrategy="uniform", condition="nTaxa", taxa=species)
print(spTree)
# let's pick a constant effective population size of 50:
popSize <- 50
# let's simulate gene trees now:
# taxa names:
for (g in 1:n_genes) {
  for (i in 1:n_species) {
    for (j in 1:n_alleles) {
        taxons[g][(i-1)*n_alleles+j] <- taxon(taxonName="Species_"+i+"_"+j, speciesName="Species_"+i)
    }
  }
  geneTrees[g] ~ dnMultiSpeciesCoalescentUniformPrior(speciesTree=spTree, max=popSize, taxa=taxons[g])
  print(geneTrees[g])
}
# We can save the species tree and the gene trees:
write(spTree, filename=dataFolder+"speciesTree")
# Saving the gene trees
for (i in 1:(n_genes)) {
  write(geneTrees[i], filename=dataFolder+"geneTree_"+i+".tree")
})");
	help_strings[string("dnMultiSpeciesCoalescentUniformPrior")][string("name")] = string(R"(dnMultiSpeciesCoalescentUniformPrior)");
	help_references[string("dnMultiSpeciesCoalescentUniformPrior")].push_back(RbHelpReference(R"(Bayes Estimation of Species Divergence Times and Ancestral Population Sizes Using DNA Sequences From Multiple Loci. Bruce Rannala and Ziheng Yang. GENETICS August 1, 2003 vol. 164 no. 4 1645-1656.)",R"()",R"(http://www.genetics.org/content/164/4/1645.short )"));
	help_references[string("dnMultiSpeciesCoalescentUniformPrior")].push_back(RbHelpReference(R"('Bayesian Inference of Species Trees from Multilocus Data. Heled and Drummond Mol. Biol Evol. 27 (3): 570-580, 2010.')",R"('DOI: https://doi.org/10.1093/molbev/msp274')",R"(https://academic.oup.com/mbe/article/27/3/570/999753/Bayesian-Inference-of-Species-Trees-from )"));
	help_references[string("dnMultiSpeciesCoalescentUniformPrior")].push_back(RbHelpReference(R"(Integration within the Felsenstein equation for improved Markov chain Monte Carlo methods in population genetics. Jody Hey and Rasmus Nielsen. PNAS. 104 (8): 2785-2790, 2007.)",R"('DOI: https://doi.org/10.1073/pnas.0611164104')",R"(https://www.pnas.org/content/104/8/2785 )"));
	help_arrays[string("dnMultiSpeciesCoalescentUniformPrior")][string("see_also")].push_back(string(R"(dnMultiSpeciesCoalescent)"));
	help_arrays[string("dnMultiSpeciesCoalescentUniformPrior")][string("see_also")].push_back(string(R"(dnMultiSpeciesCoalescentInverseGamma)"));
	help_strings[string("dnMultiSpeciesCoalescentUniformPrior")][string("title")] = string(R"(Multispecies coalescent Distribution)");
	help_arrays[string("dnMultiValueEvent")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("dnMultiValueEvent")][string("description")] = string(R"(A MultiValueEvent distribution.)");
	help_strings[string("dnMultiValueEvent")][string("name")] = string(R"(dnMultiValueEvent)");
	help_arrays[string("dnMultiValueEvent")][string("see_also")].push_back(string(R"(dnEvent)"));
	help_strings[string("dnMultiValueEvent")][string("title")] = string(R"(MultiValueEvent Distribution)");
	help_arrays[string("dnMultinomial")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("dnMultinomial")][string("description")] = string(R"(A multinomial distribution defines a probability distribution on a vector of natural numbers. It is understood as randomly picking n times from the k categories with replacement where each catefory has its own probability p[i].)");
	help_strings[string("dnMultinomial")][string("example")] = string(R"(p <- simplex(1,1,1,1)
x ~ dnMultinomial(10, p)
y ~ dnDirichlet(x)
y.clamp( simplex(1,2,3,4) )
moves[1] = mvSlide(x, delta=0.1, weight=1.0)
monitors[1] = screenmonitor(printgen=1000, separator = "        ", x)
mymodel = model(p)
mymcmc = mcmc(mymodel, monitors, moves)
mymcmc.burnin(generations=20000,tuningInterval=100)
mymcmc.run(generations=200000))");
	help_strings[string("dnMultinomial")][string("name")] = string(R"(dnMultinomial)");
	help_arrays[string("dnMultinomial")][string("see_also")].push_back(string(R"(dnDirichlet)"));
	help_strings[string("dnMultinomial")][string("title")] = string(R"(Multinomial Distribution)");
	help_arrays[string("dnMultivariateNormal")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_arrays[string("dnMultivariateNormal")][string("authors")].push_back(string(R"(Michael Landis)"));
	help_strings[string("dnMultivariateNormal")][string("description")] = string(R"(The multivariate normal distribution has the probability density:

f(x) = det(2 pi Sigma)^(-1/2) e^{-(1/2) (x-mu)' Sigma^-1 (x-mu)}

where mu is a vector of mean values and Sigma is a covariance matrix. Note, this distribution may also be parameterized in terms of the precision matrix, Sigma^-1.)");
	help_strings[string("dnMultivariateNormal")][string("example")] = string(R"(dim = 4
df = 100
kappa <- 2
Sigma ~ dnWishart(df, kappa, dim)
for (i in 1:dim) { mu[i] ~ dnUnif(-1, 1) }
x ~ dnMultivariateNormal( mean=mu, covariance=Sigma )
mv[1] = mvCorrelationMatrixElementSwap(Sigma)
mv[2] = mvCorrelationMatrixRandomWalk(Sigma)
mv[3] = mvCorrelationMatrixSingleElementBeta(Sigma)
mv[4] = mvCorrelationMatrixSpecificElementBeta(Sigma)
mv[5] = mvCorrelationMatrixUpdate(Sigma)
mv[6] = mvVectorSlide(x))");
	help_strings[string("dnMultivariateNormal")][string("name")] = string(R"(dnMultivariateNormal)");
	help_arrays[string("dnMultivariateNormal")][string("see_also")].push_back(string(R"(dnNormal)"));
	help_arrays[string("dnMultivariateNormal")][string("see_also")].push_back(string(R"(dnWishart)"));
	help_arrays[string("dnMultivariateNormal")][string("see_also")].push_back(string(R"(mvCorrelationMatrixUpdate)"));
	help_strings[string("dnMultivariateNormal")][string("title")] = string(R"(Multivariate Normal Distribution)");
	help_arrays[string("dnNbinomial")][string("authors")].push_back(string(R"(Walker Pett)"));
	help_strings[string("dnNbinomial")][string("description")] = string(R"(Negative binomial probability distribution of x successes before r failures.)");
	help_strings[string("dnNbinomial")][string("details")] = string(R"(The negative binomial probability distribution describes the number of successes before r failures, where the success probability is p. The probability is given by (x + r - 1 choose x) p^(x) * (1-p)^r)");
	help_strings[string("dnNbinomial")][string("example")] = string(R"(p ~ dnBeta(1.0,1.0)
x ~ dnNegativeBinomial(r=10,p)
x.clamp(8)
moves[1] = mvSlide(p, delta=0.1, weight=1.0)
monitors[1] = screenmonitor(printgen=1000, separator = "        ", x)
mymodel = model(p)
mymcmc = mcmc(mymodel, monitors, moves)
mymcmc.burnin(generations=20000,tuningInterval=100)
mymcmc.run(generations=200000))");
	help_strings[string("dnNbinomial")][string("name")] = string(R"(dnNbinomial)");
	help_arrays[string("dnNbinomial")][string("see_also")].push_back(string(R"(dnBinomial)"));
	help_strings[string("dnNbinomial")][string("title")] = string(R"(Negative Binomial Distribution)");
	help_arrays[string("dnNormal")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("dnNormal")][string("description")] = string(R"(Normal (gaussian) distribution with mean equal to ‘mean’ and standard deviation equal to ‘sd’.)");
	help_strings[string("dnNormal")][string("details")] = string(R"(The normal distribution has density:

f(x) = 1/(sqrt(2 pi) sigma) e^-((x - mu)^2/(2 sigma^2))

where mu is the mean of the distribution and sigma the standard deviation.)");
	help_strings[string("dnNormal")][string("example")] = string(R"(# we simulate some observations
x <- rnorm(n=10,mean=5,sd=10)
# let's see what the minimum is (you could do the max too)
min(x)
# let's also see what the mean and the variance are
mean(x)
var(x)
sd(x))");
	help_strings[string("dnNormal")][string("name")] = string(R"(dnNormal)");
	help_arrays[string("dnNormal")][string("see_also")].push_back(string(R"(dnLognormal)"));
	help_strings[string("dnNormal")][string("title")] = string(R"(Normal Distribution)");
	help_arrays[string("dnOrnsteinUhlenbeck")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("dnOrnsteinUhlenbeck")][string("description")] = string(R"(A Bernoulli-distributed random variable takes the value 1 with probability p and the value 0 with probability 1-p.)");
	help_strings[string("dnOrnsteinUhlenbeck")][string("example")] = string(R"(p ~ dnBeta(1.0,1.0)
x ~ dnBernoulli(p)
x.clamp(1)
moves[1] = mvSlide(p, delta=0.1, weight=1.0)
monitors[1] = screenmonitor(printgen=1000, separator = "        ", speciation)
mymodel = model(p)
mymcmc = mcmc(mymodel, monitors, moves)
mymcmc.burnin(generations=20000,tuningInterval=100)
mymcmc.run(generations=200000))");
	help_strings[string("dnOrnsteinUhlenbeck")][string("name")] = string(R"(dnOrnsteinUhlenbeck)");
	help_arrays[string("dnOrnsteinUhlenbeck")][string("see_also")].push_back(string(R"(dnBinomial)"));
	help_strings[string("dnOrnsteinUhlenbeck")][string("title")] = string(R"(Bernoulli Distribution)");
	help_strings[string("dnOutgroupBirthDeath")][string("name")] = string(R"(dnOutgroupBirthDeath)");
	help_strings[string("dnPhyloBrownian")][string("name")] = string(R"(dnPhyloBrownian)");
	help_strings[string("dnPhyloBrownianMVN")][string("name")] = string(R"(dnPhyloBrownianMVN)");
	help_strings[string("dnPhyloBrownianMultiSampleREML")][string("name")] = string(R"(dnPhyloBrownianMultiSampleREML)");
	help_strings[string("dnPhyloBrownianMultiVariate")][string("name")] = string(R"(dnPhyloBrownianMultiVariate)");
	help_strings[string("dnPhyloBrownianREML")][string("name")] = string(R"(dnPhyloBrownianREML)");
	help_strings[string("dnPhyloCTMC")][string("name")] = string(R"(dnPhyloCTMC)");
	help_strings[string("dnPhyloCTMCClado")][string("name")] = string(R"(dnPhyloCTMCClado)");
	help_strings[string("dnPhyloCTMCDASequence")][string("name")] = string(R"(dnPhyloCTMCDASequence)");
	help_strings[string("dnPhyloCTMCDASiteIID")][string("name")] = string(R"(dnPhyloCTMCDASiteIID)");
	help_strings[string("dnPhyloCTMCDollo")][string("name")] = string(R"(dnPhyloCTMCDollo)");
	help_strings[string("dnPhyloDistanceGamma")][string("name")] = string(R"(dnPhyloDistanceGamma)");
	help_strings[string("dnPhyloMultiSampleOrnsteinUhlenbeck")][string("name")] = string(R"(dnPhyloMultiSampleOrnsteinUhlenbeck)");
	help_strings[string("dnPhyloMultiSampleOrnsteinUhlenbeckREML")][string("name")] = string(R"(dnPhyloMultiSampleOrnsteinUhlenbeckREML)");
	help_strings[string("dnPhyloMultivariateBrownianMultiSampleREML")][string("name")] = string(R"(dnPhyloMultivariateBrownianMultiSampleREML)");
	help_arrays[string("dnPhyloMultivariateBrownianREML")][string("authors")].push_back(string(R"(Michael R. May)"));
	help_arrays[string("dnPhyloMultivariateBrownianREML")][string("authors")].push_back(string(R"(Nicolai Vetr)"));
	help_strings[string("dnPhyloMultivariateBrownianREML")][string("description")] = string(R"(Multivariate Brownian motion over a phylogeny with variance-covariance matrix rateMatrix.)");
	help_strings[string("dnPhyloMultivariateBrownianREML")][string("example")] = string(R"(
# generate a tree and variance-covariance matrix.
psi ~ dnUniformTimeTree(1, [taxon("A"),taxon("B"),taxon("C")])
Sigma <- diagonalMatrix(5)

# generate the multivariate data.
x ~ dnPhyloMultivariateBrownianREML(tree=psi, rateMatrix=Sigma)

# print the simulated data.
x)");
	help_strings[string("dnPhyloMultivariateBrownianREML")][string("name")] = string(R"(dnPhyloMultivariateBrownianREML)");
	help_references[string("dnPhyloMultivariateBrownianREML")].push_back(RbHelpReference(R"(Huelsenbeck JP, B Rannala (2003). Detecting correlation between characters in a comparative analysis with uncertain phylogeny. Evolution, 57(6):1237-1247.)",R"()",R"()"));
	help_arrays[string("dnPhyloMultivariateBrownianREML")][string("see_also")].push_back(string(R"(dnPhyloBrownianREML)"));
	help_arrays[string("dnPhyloMultivariateBrownianREML")][string("see_also")].push_back(string(R"(dnPhyloBrownianMVN)"));
	help_strings[string("dnPhyloMultivariateBrownianREML")][string("title")] = string(R"(Phylogenetic Multivariate Brownian Motion)");
	help_strings[string("dnPhyloOrnsteinUhlenbeck")][string("name")] = string(R"(dnPhyloOrnsteinUhlenbeck)");
	help_strings[string("dnPhyloOrnsteinUhlenbeckMVN")][string("name")] = string(R"(dnPhyloOrnsteinUhlenbeckMVN)");
	help_strings[string("dnPhyloOrnsteinUhlenbeckREML")][string("name")] = string(R"(dnPhyloOrnsteinUhlenbeckREML)");
	help_strings[string("dnPhyloOrnsteinUhlenbeckThreePoint")][string("name")] = string(R"(dnPhyloOrnsteinUhlenbeckThreePoint)");
	help_strings[string("dnPhyloWhiteNoise")][string("name")] = string(R"(dnPhyloWhiteNoise)");
	help_arrays[string("dnPointMass")][string("authors")].push_back(string(R"(Walker Pett)"));
	help_strings[string("dnPointMass")][string("description")] = string(R"(Point mass distribution.)");
	help_strings[string("dnPointMass")][string("details")] = string(R"(The point mass distribution, or Dirac delta function, has density f(x) = 1 when x is equal to the point mass value.)");
	help_strings[string("dnPointMass")][string("example")] = string(R"(u ~ dnPointMass(1.2))");
	help_strings[string("dnPointMass")][string("name")] = string(R"(dnPointMass)");
	help_strings[string("dnPointMass")][string("title")] = string(R"(Point Mass Distribution)");
	help_arrays[string("dnPoisson")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("dnPoisson")][string("description")] = string(R"(A Poisson distribution defines probabilities for natural numbers. It is defined as the number of exponentially distributed events in a given interval.)");
	help_strings[string("dnPoisson")][string("example")] = string(R"(l ~ dnUniform(0.0,100.0)
x ~ dnPoisson(l)
x.clamp(10)
moves[1] = mvSlide(l, delta=0.1, weight=1.0)
monitors[1] = mnScreen(printgen=1000, separator = "        ", l)
mymodel = model(l)
mymcmc = mcmc(mymodel, monitors, moves)
mymcmc.burnin(generations=20000,tuningInterval=100)
mymcmc.run(generations=200000))");
	help_strings[string("dnPoisson")][string("name")] = string(R"(dnPoisson)");
	help_arrays[string("dnPoisson")][string("see_also")].push_back(string(R"(dnGeom)"));
	help_strings[string("dnPoisson")][string("title")] = string(R"(Poisson Distribution)");
	help_strings[string("dnReversibleJumpMixture")][string("name")] = string(R"(dnReversibleJumpMixture)");
	help_strings[string("dnSBBDP")][string("name")] = string(R"(dnSBBDP)");
	help_arrays[string("dnScaledDirichlet")][string("authors")].push_back(string(R"(Andrew Magee)"));
	help_strings[string("dnScaledDirichlet")][string("description")] = string(R"(Scaled Dirichlet probability distribution on a simplex.)");
	help_strings[string("dnScaledDirichlet")][string("details")] = string(R"(The scaled Dirichlet probability distribution is the generalization of the dirichlet distribution. A random variable from a scaled Dirichlet distribution is a simplex, i.e., a vector of probabilities that sum to 1. If b[1]=b[2]=...=b[n], then the scaledDirichlet(alpha,beta) collapses to the Dirichlet with the same alphas.)");
	help_strings[string("dnScaledDirichlet")][string("example")] = string(R"(# lets get a draw from a Dirichlet distribution
a <- [1,1,1,1]   # we could also use rep(1,4)
b <- [1,2,3,4]   # if these are all equal, the scaled Dirichlet is equivilent to the Dirichlet(a)x ~ dnScaledDirichlet(a,b)
x
# let check if b really sums to 1
sum(x))");
	help_strings[string("dnScaledDirichlet")][string("name")] = string(R"(dnScaledDirichlet)");
	help_arrays[string("dnScaledDirichlet")][string("see_also")].push_back(string(R"(dnDirichlet)"));
	help_arrays[string("dnScaledDirichlet")][string("see_also")].push_back(string(R"(simplex)"));
	help_strings[string("dnScaledDirichlet")][string("title")] = string(R"(Scaled Dirichlet Distribution)");
	help_strings[string("dnSerialSampledBirthDeath")][string("name")] = string(R"(dnSerialSampledBirthDeath)");
	help_arrays[string("dnSoftBoundUniformNormal")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("dnSoftBoundUniformNormal")][string("description")] = string(R"(A softbound uniform distribution with normally distributed tails outside the interval of the uniform distribution.)");
	help_strings[string("dnSoftBoundUniformNormal")][string("details")] = string(R"(The center piece of this distribution a uniform distribution defined between the given interval. A variable is drawn from that uniform distribution with probability p and with probability 1-p outside the interval. The probability density outside the interval is described by a normal distribution with standard deviation sd.)");
	help_strings[string("dnSoftBoundUniformNormal")][string("example")] = string(R"(p ~ dnBeta(1.0,1.0)
x ~ dnBernoulli(p)
x.clamp(1)
moves[1] = mvSlide(p, delta=0.1, weight=1.0)
monitors[1] = screenmonitor(printgen=1000, separator = "        ", speciation)
mymodel = model(p)
mymcmc = mcmc(mymodel, monitors, moves)
mymcmc.burnin(generations=20000,tuningInterval=100)
mymcmc.run(generations=200000))");
	help_strings[string("dnSoftBoundUniformNormal")][string("name")] = string(R"(dnSoftBoundUniformNormal)");
	help_arrays[string("dnSoftBoundUniformNormal")][string("see_also")].push_back(string(R"(dnUniform)"));
	help_strings[string("dnSoftBoundUniformNormal")][string("title")] = string(R"(Softbound Uniform Distribution with Normal distributed tails.)");
	help_arrays[string("dnStudentT")][string("authors")].push_back(string(R"(Wade Dismukes and Kevin Quinteros)"));
	help_strings[string("dnStudentT")][string("description")] = string(R"(The student's t probability distribution.)");
	help_strings[string("dnStudentT")][string("example")] = string(R"(# The most important use of the Student T distribution
# is arguable the quantile function.
# You can access it the following way:
df <- 10
a := qStudentT(0.025, df)
a)");
	help_strings[string("dnStudentT")][string("name")] = string(R"(dnStudentT)");
	help_strings[string("dnStudentT")][string("title")] = string(R"(Student T Distribution)");
	help_strings[string("dnTimeVaryingStateDependentSpeciationExtinction")][string("name")] = string(R"(dnTimeVaryingStateDependentSpeciationExtinction)");
	help_strings[string("dnUPP")][string("name")] = string(R"(dnUPP)");
	help_strings[string("dnUltrametricTree")][string("name")] = string(R"(dnUltrametricTree)");
	help_arrays[string("dnUniform")][string("authors")].push_back(string(R"(Michael Landis)"));
	help_strings[string("dnUniform")][string("description")] = string(R"(Uniform distribution with lower and uppper bounds.)");
	help_strings[string("dnUniform")][string("details")] = string(R"(The uniform distribution has density, f(x) = 1/(b-a), where b is the upper bound and a is the lower bound.)");
	help_strings[string("dnUniform")][string("example")] = string(R"(# set the lower bound
a <- -2.5
# set the upper bound
b <- -3.9
# create a stochastic node with a uniform prior
u ~ dnUniform(a, b))");
	help_strings[string("dnUniform")][string("name")] = string(R"(dnUniform)");
	help_strings[string("dnUniform")][string("title")] = string(R"(Uniform Distribution)");
	help_strings[string("dnUniformInteger")][string("name")] = string(R"(dnUniformInteger)");
	help_strings[string("dnUniformNatural")][string("name")] = string(R"(dnUniformNatural)");
	help_strings[string("dnUniformTimeTree")][string("name")] = string(R"(dnUniformTimeTree)");
	help_strings[string("dnUniformTopology")][string("name")] = string(R"(dnUniformTopology)");
	help_strings[string("dnUniformTopologyBranchLength")][string("name")] = string(R"(dnUniformTopologyBranchLength)");
	help_arrays[string("dnVarianceGamma")][string("authors")].push_back(string(R"(Michael Landis)"));
	help_strings[string("dnVarianceGamma")][string("description")] = string(R"(Variance-gamma distribution with location ‘mu’.)");
	help_strings[string("dnVarianceGamma")][string("details")] = string(R"(The variance-gamma distribution has density:

f(x) = 1/(sqrt(2 pi) sigma) e^-((x - mu)^2/(2 sigma^2))

where mu is the mean of the distribution and sigma the standard deviation.)");
	help_strings[string("dnVarianceGamma")][string("example")] = string(R"(# we simulate some oversations
x <- rnorm(n=10,mean=5,sd=10)
# let's see what the minum is (you could do the max too)
min(x)
# let's also see what the mean and the variance are
mean(x)
var(x)
sd(x))");
	help_strings[string("dnVarianceGamma")][string("name")] = string(R"(dnVarianceGamma)");
	help_strings[string("dnVarianceGamma")][string("title")] = string(R"(Variance-gamma Distribution)");
	help_strings[string("dnWeightedConstrainedNodeOrder")][string("name")] = string(R"(dnWeightedConstrainedNodeOrder)");
	help_strings[string("dnWeightedSample")][string("name")] = string(R"(dnWeightedSample)");
	help_arrays[string("dnWhiteNoise")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("dnWhiteNoise")][string("description")] = string(R"(White-Noise process for positive real numbers.)");
	help_strings[string("dnWhiteNoise")][string("details")] = string(R"(The white-noise process is a process of a positive continuous variable similar to Brownian motion and the Ornstein-Uhlenbeck process. However, the white-noise process has a large variance when the time is small, and has small variance if the time is large.)");
	help_strings[string("dnWhiteNoise")][string("example")] = string(R"(# lets simulate
a <- rWhiteNoise(1000,mu=1,sigma=4,time=4)
# we expect a mean of 1
mean(a)

# create a random variable
x ~ dnWhiteNoise(mu=1.0,sigma=4,time=1)
x)");
	help_strings[string("dnWhiteNoise")][string("name")] = string(R"(dnWhiteNoise)");
	help_arrays[string("dnWhiteNoise")][string("see_also")].push_back(string(R"(dnGamma)"));
	help_strings[string("dnWhiteNoise")][string("title")] = string(R"(White-Noise Process)");
	help_strings[string("dnWishart")][string("name")] = string(R"(dnWishart)");
	help_arrays[string("exists")][string("authors")].push_back(string(R"(Michael Landis)"));
	help_strings[string("exists")][string("description")] = string(R"(Determines whether the RevBayes workspace contains a variable named 'name')");
	help_strings[string("exists")][string("details")] = string(R"('exists' returns 'true' if the workspace contains a variable whose name matches the String 'name' and 'false' otherwise. One use of 'exists' is to add Move and Monitor objects conditional on the variable 'x' existing. The function 'ls' provides a summary for all variable names that 'exists' would evaluate as 'true'.)");
	help_strings[string("exists")][string("example")] = string(R"(## Correct usage: does "x" exist?
x <- 1.0
exists("x")

## Incorrect usage: does "1.0" exist?
exists(x))");
	help_strings[string("exists")][string("name")] = string(R"(exists)");
	help_arrays[string("exists")][string("see_also")].push_back(string(R"(clear)"));
	help_strings[string("exists")][string("title")] = string(R"(Check whether a variable exists)");
	help_strings[string("exp")][string("name")] = string(R"(exp)");
	help_strings[string("floor")][string("name")] = string(R"(floor)");
	help_strings[string("fnAdjacentRateModifier")][string("name")] = string(R"(fnAdjacentRateModifier)");
	help_strings[string("fnBetaBrokenStick")][string("name")] = string(R"(fnBetaBrokenStick)");
	help_strings[string("fnBinaryMutationCoalescentRateMatrix")][string("name")] = string(R"(fnBinaryMutationCoalescentRateMatrix)");
	help_strings[string("fnBiogeoDE")][string("name")] = string(R"(fnBiogeoDE)");
	help_strings[string("fnBiogeographyCladoEventsBD")][string("name")] = string(R"(fnBiogeographyCladoEventsBD)");
	help_strings[string("fnBlosum62")][string("name")] = string(R"(fnBlosum62)");
	help_strings[string("fnChromosomes")][string("name")] = string(R"(fnChromosomes)");
	help_strings[string("fnChromosomesCladoEventsBD")][string("name")] = string(R"(fnChromosomesCladoEventsBD)");
	help_strings[string("fnChromosomesCladoProbs")][string("name")] = string(R"(fnChromosomesCladoProbs)");
	help_strings[string("fnChromosomesPloidy")][string("name")] = string(R"(fnChromosomesPloidy)");
	help_strings[string("fnChromosomesPloidyCladoEventsBD")][string("name")] = string(R"(fnChromosomesPloidyCladoEventsBD)");
	help_strings[string("fnChronoToPhylo")][string("name")] = string(R"(fnChronoToPhylo)");
	help_strings[string("fnCladeSpecificHierarchicalBranchRate")][string("name")] = string(R"(fnCladeSpecificHierarchicalBranchRate)");
	help_strings[string("fnCladogeneticSpeciationRateMatrix")][string("name")] = string(R"(fnCladogeneticSpeciationRateMatrix)");
	help_arrays[string("fnCoala")][string("authors")].push_back(string(R"(Bastien Boussau)"));
	help_strings[string("fnCoala")][string("name")] = string(R"(fnCoala)");
	help_references[string("fnCoala")].push_back(RbHelpReference(R"(A branch-heterogeneous model of protein evolution for efficient inference of ancestral sequences. Groussin M, Boussau B, Gouy M. Syst Biol. 2013 Jul;62(4):523-38.)",R"(10.1093/sysbio/syt016)",R"(https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3676677/ )"));
	help_strings[string("fnCodon")][string("name")] = string(R"(fnCodon)");
	help_strings[string("fnCodonGY94")][string("description")] = string(R"(The Goldman-Yang (1994) codon model.

A rate matrix on the 61 non-stop codons (in the standard genetic code).

Rates between codons with more than one nucleotide change are equal to 0.

In this model the rate Q(i,j) from i -> j is proportional to the frequency of codon j.
This means that the rate of change between low-frequency codons is lower than the rate
between high-frequency codons, even when the nucleotide change involved is the same.
In other words, the rate of change from nucleotide n1 -> n2 depends on its neighboring
nucleotides.  This differs from the Muse-Gaut (1994) model, and is perhaps less realistic.

Unlike the Muse-Gaut (1994) model, the Goldman-Yang (1994) model can allow all the codon
frequencies to vary independently.)");
	help_strings[string("fnCodonGY94")][string("example")] = string(R"(kappa ~ dnLognormal(0,1)
omega ~ dnUniform(0,1)
pi61 ~ dnDirichlet( rep(2.0, 61) )
Q1 := fnCodonGY94( kappa, omega, pi61 )

pi1 ~ dnDirichlet( rep(2.0, 4) )
Q2 := fnCodonGY94( kappa, omega, fnF1x4(pi1) )

pi2 ~ dnDirichlet( rep(2.0, 4) )
pi3 ~ dnDirichlet( rep(2.0, 4) )
Q3 := fnCodonGY94( kappa, omega, fnF3x4(pi1, pi2, pi3) ))");
	help_strings[string("fnCodonGY94")][string("name")] = string(R"(fnCodonGY94)");
	help_references[string("fnCodonGY94")].push_back(RbHelpReference(R"(Goldman, N. and Z. Yang (1994). A codon-based model of nucleotide substitution for protein-coding DNA sequences. Mol. Biol. Evol. (1994) 11 (5):725-736)",R"(https://doi.org/10.1093/oxfordjournals.molbev.a040153 )",R"()"));
	help_arrays[string("fnCodonGY94")][string("see_also")].push_back(string(R"(fnF1x4, fnF3x4, fnCodonMG94, fnCodonMG94K)"));
	help_strings[string("fnCodonGY94")][string("title")] = string(R"(The Goldman-Yang (1994) codon rate matrix)");
	help_strings[string("fnCodonHKY")][string("name")] = string(R"(fnCodonHKY)");
	help_strings[string("fnCodonMG94")][string("description")] = string(R"(The Muse-Gaut (1994) codon model.

A rate matrix on the 61 non-stop codons (in the standard genetic code).

Rates between codons with more than one nucleotide change are equal to 0.

In this model the rate Q(i,j) from i -> j is proportional to the frequency of
nucleotide in codon j that changed.  This differs from the Goldman-Yang (1994) model,
where Q(i,j) is proportional to the frequency of the entire codon j.

Unlike the Goldman-Yang (1994) model, the Muse-Gaut (1994) model does not allow all the codon
frequencies to vary independently.)");
	help_strings[string("fnCodonMG94")][string("example")] = string(R"(omega ~ dnUniform(0,1)
pi ~ dnDirichlet( rep(2.0, 4) )
Q1 := fnCodonMG94( omega, pi )

Q2 := fndNdS( omega, fnX3( fnF81(pi) ) ) # MG94 = F81 + X3 + dNdS)");
	help_strings[string("fnCodonMG94")][string("name")] = string(R"(fnCodonMG94)");
	help_references[string("fnCodonMG94")].push_back(RbHelpReference(R"(Muse, S. and B. Gaut (1994) A likelihood approach for comparing synonymous and nonsynonymous nucleotide substitution rates, with application to the chloroplast genome. Mol. Biol. Evol. (1994) 11 (5):715-724)",R"(https://doi.org/10.1093/oxfordjournals.molbev.a040152 )",R"()"));
	help_arrays[string("fnCodonMG94")][string("see_also")].push_back(string(R"(fnCodonMG94, fnCodonMG94K)"));
	help_strings[string("fnCodonMG94")][string("title")] = string(R"(The Muse-Gaut (1994) codon rate matrix)");
	help_strings[string("fnCodonMG94K")][string("description")] = string(R"(The Muse-Gaut (1994) codon model, extended with a transition/transversion rate ratio.

A rate matrix on the 61 non-stop codons (in the standard genetic code).

Rates between codons with more than one nucleotide change are equal to 0.

In this model the rate Q(i,j) from i -> j is proportional to the frequency of
nucleotide in codon j that changed.  This differs from the Goldman-Yang (1994) model,
where Q(i,j) is proportional to the frequency of the entire codon j.

This version is an extension of the fnCodonMG94 model to add a transition/transversion
rate ratio.  This makes it more comparable to the Goldman-Yang (1994) model.

Unlike the Goldman-Yang (1994) model, the Muse-Gaut (1994) model does not allow all the codon
frequencies to vary independently.)");
	help_strings[string("fnCodonMG94K")][string("example")] = string(R"(kappa ~ dnLognormal(0,1)
omega ~ dnUniform(0,1)
pi ~ dnDirichlet( rep(2.0, 4) )
Q1 := fnCodonMG94K( kappa, omega, pi )

Q2 := fndNdS( omega, fnX3( fnHKY( kappa, pi) ) ) # MG94K = HKY + X3 + dNdS)");
	help_strings[string("fnCodonMG94K")][string("name")] = string(R"(fnCodonMG94K)");
	help_references[string("fnCodonMG94K")].push_back(RbHelpReference(R"(Muse, S. and B. Gaut (1994) A likelihood approach for comparing synonymous and nonsynonymous nucleotide substitution rates, with application to the chloroplast genome. Mol. Biol. Evol. (1994) 11 (5):715-724)",R"(https://doi.org/10.1093/oxfordjournals.molbev.a040152 )",R"()"));
	help_arrays[string("fnCodonMG94K")][string("see_also")].push_back(string(R"(fnCodonGY94, fnCodonMG94K)"));
	help_strings[string("fnCodonMG94K")][string("title")] = string(R"(The Muse-Gaut (1994) codon rate matrix + K.)");
	help_strings[string("fnCovarion")][string("name")] = string(R"(fnCovarion)");
	help_strings[string("fnCovarionRateMatrix")][string("name")] = string(R"(fnCovarionRateMatrix)");
	help_strings[string("fnCpRev")][string("name")] = string(R"(fnCpRev)");
	help_strings[string("fnDECCladoProbs")][string("name")] = string(R"(fnDECCladoProbs)");
	help_strings[string("fnDECRateMatrix")][string("name")] = string(R"(fnDECRateMatrix)");
	help_strings[string("fnDECRates")][string("name")] = string(R"(fnDECRates)");
	help_strings[string("fnDECRoot")][string("name")] = string(R"(fnDECRoot)");
	help_strings[string("fnDayhoff")][string("name")] = string(R"(fnDayhoff)");
	help_strings[string("fnDecompVarCovar")][string("name")] = string(R"(fnDecompVarCovar)");
	help_strings[string("fnDiscretizeBeta")][string("name")] = string(R"(fnDiscretizeBeta)");
	help_strings[string("fnDiscretizeBetaQuadrature")][string("name")] = string(R"(fnDiscretizeBetaQuadrature)");
	help_strings[string("fnDiscretizeDistribution")][string("name")] = string(R"(fnDiscretizeDistribution)");
	help_strings[string("fnDiscretizeGamma")][string("name")] = string(R"(fnDiscretizeGamma)");
	help_strings[string("fnDiscretizeGammaFromBetaQuantiles")][string("name")] = string(R"(fnDiscretizeGammaFromBetaQuantiles)");
	help_strings[string("fnDiscretizeGammaQuadrature")][string("name")] = string(R"(fnDiscretizeGammaQuadrature)");
	help_strings[string("fnDiscretizeLognormalQuadrature")][string("name")] = string(R"(fnDiscretizeLognormalQuadrature)");
	help_strings[string("fnDistanceRateModifier")][string("name")] = string(R"(fnDistanceRateModifier)");
	help_strings[string("fnDppConcFromMean")][string("name")] = string(R"(fnDppConcFromMean)");
	help_strings[string("fnDppMeanFromConc")][string("name")] = string(R"(fnDppMeanFromConc)");
	help_strings[string("fnEarlyBurst")][string("name")] = string(R"(fnEarlyBurst)");
	help_strings[string("fnEpoch")][string("name")] = string(R"(fnEpoch)");
	help_strings[string("fnEpochCladoProbs")][string("name")] = string(R"(fnEpochCladoProbs)");
	help_strings[string("fnExtantTree")][string("name")] = string(R"(fnExtantTree)");
	help_strings[string("fnF1x4")][string("description")] = string(R"(This treats codon frequencies as a product of independent nucleotide frequencies.

Since stop codons are removed from the codon alphabet, frequencies are renormalized
so that the frequencies of non-stop codons sum to 1.0.)");
	help_strings[string("fnF1x4")][string("example")] = string(R"(kappa ~ dnLognormal(0,1)
omega ~ dnUniform(0,1)
pi ~ dnDirichlet( v(2.0, 2.0, 2.0, 2.0) )
Q := fnCodonGY94( kappa, omega, fnF1x4(pi) ))");
	help_strings[string("fnF1x4")][string("name")] = string(R"(fnF1x4)");
	help_arrays[string("fnF1x4")][string("see_also")].push_back(string(R"(fnGY94, fnF3x4)"));
	help_strings[string("fnF1x4")][string("title")] = string(R"(The F1x4 codon frequency model)");
	help_strings[string("fnF2x4")][string("description")] = string(R"(This treats doublet frequencies as a product of independent nucleotide frequencies.)");
	help_strings[string("fnF2x4")][string("example")] = string(R"(# An RNA stem model with independent base frequencies (from fnF2x4),
# and simultaneous 2-nucleotide changes allows.
nuc_pi ~ dnDirichlet( v(2.0, 2.0, 2.0, 2.0) )
rna_stem_er ~ dnDirichlet( rep(1.0, 16*15/2) )
rna_stem_pi := fnF2x4(nuc_pi, nuc_pi)
Q := fnGTR(rna_stem_er, rna_stem_pi))");
	help_strings[string("fnF2x4")][string("name")] = string(R"(fnF2x4)");
	help_arrays[string("fnF2x4")][string("see_also")].push_back(string(R"(fnX2)"));
	help_strings[string("fnF2x4")][string("title")] = string(R"(The F2x4 doublet frequency model)");
	help_strings[string("fnF3x4")][string("description")] = string(R"(This treats codon frequencies as a product of independent nucleotide frequencies.

Since stop codons are removed from the codon alphabet, frequencies are renormalized
so that the frequencies of non-stop codons sum to 1.0.)");
	help_strings[string("fnF3x4")][string("example")] = string(R"(kappa ~ dnLognormal(0,1)
omega ~ dnUniform(0,1)
pi1 ~ dnDirichlet( v(2.0, 2.0, 2.0, 2.0) )
pi2 ~ dnDirichlet( v(2.0, 2.0, 2.0, 2.0) )
pi3 ~ dnDirichlet( v(2.0, 2.0, 2.0, 2.0) )
Q := fnCodonGY94( kappa, omega, fnF3x4(pi1, pi2, pi3) ))");
	help_strings[string("fnF3x4")][string("name")] = string(R"(fnF3x4)");
	help_arrays[string("fnF3x4")][string("see_also")].push_back(string(R"(fnGY94, fnF1x4)"));
	help_strings[string("fnF3x4")][string("title")] = string(R"(The F3x4 codon frequency model)");
	help_strings[string("fnF81")][string("name")] = string(R"(fnF81)");
	help_strings[string("fnFMutSel")][string("description")] = string(R"(Constructs a rate matrix from 61 scaled selection coefficients w[i] and
a 4x4 nucleotide mutation rate matrix mu(i,j).  In the original paper
the nucleotide mutation rate matrix is a GTR rate matrix.

The FMutSel0 model differs from FMutSel by constraining all codons for
the same amino acid to have the same scaled selection coefficient.

The function fnMutSel differs from fnFMutSel by taking a codon mutation
rate matrix.

A substitution from allele i -> j can be decomposed into
 (1) all individuals initially have state i
 (2) a single individual mutates from i -> j, at rate mu(i,j)
 (3) the allele j goes to fixation

Then the substitution rate Q is then given by
  Q(i,j) = mu(i,j) * Pr(j goes to fixation | i was fixed previously).

The probability of fixation is determined by scaled selection coefficients:
  F[i] = 2*N*s[i]
and the initial frequency 1/N of allele j.)");
	help_strings[string("fnFMutSel")][string("example")] = string(R"(er ~ dnDirichlet( v(1,1,1,1,1,1) )
nuc_pi ~ dnDirichlet( rep(2.0, 4) )
F ~ dnIID(61, dnNormal(0,1))
omega ~ dnUniform(0,1)
# The FMutSel model from Yang and Nielsen (2008)
Q1 := fnFMutSel(fnGTR(er, nuc_pi), F, omega)

# The same -- fMutSel = GTR(er,nuc_pi) + X3 + MutSel(F) + dNdS(omega)
Q2 := fndNdS(fnMutSel(F, fnX3(fnGTR(er, nuc_pi))), omega))");
	help_strings[string("fnFMutSel")][string("name")] = string(R"(fnFMutSel)");
	help_references[string("fnFMutSel")].push_back(RbHelpReference(R"(Yang, Z. and R. Nielsen. Mutation-Selection Models of Codon Substitution and Their Use to Estimate Selective Strengths on Codon Usage.  Mol. Biol. Evol. (2008) 25(3):568--579)",R"(https://doi.org/10.1093/molbev/msm284 )",R"()"));
	help_arrays[string("fnFMutSel")][string("see_also")].push_back(string(R"(fnCodonGY94, fnCodonMG94, fnFMutSel0, fnMutSel)"));
	help_strings[string("fnFMutSel")][string("title")] = string(R"(The FMutSel model)");
	help_strings[string("fnFMutSel0")][string("description")] = string(R"(Constructs a rate matrix from 61 scaled selection coefficients w[i] and
a 4x4 nucleotide mutation rate matrix mu(i,j).  In the original paper
the nucleotide mutation rate matrix is a GTR rate matrix.

The FMutSel0 model is a restriction of the FMutSel model that constrains
all codons for the same amino acid to have the same scaled selection
coefficient.

The function fnMutSelAA differs from fnFMutSel0 by taking a codon mutation
rate matrix.

A substitution from allele i -> j can be decomposed into
 (1) all individuals initially have state i
 (2) a single individual mutates from i -> j, at rate mu(i,j)
 (3) the allele j goes to fixation

Then the substitution rate Q is then given by
  Q(i,j) = mu(i,j) * Pr(j goes to fixation | i was fixed previously).

The probability of fixation is determined by scaled selection coefficients:
  F[i] = 2*N*s[i]
and the initial frequency 1/N of allele j.)");
	help_strings[string("fnFMutSel0")][string("example")] = string(R"(er ~ dnDirichlet( v(1,1,1,1,1,1) )
nuc_pi ~ dnDirichlet( rep(2.0, 4) )
F ~ dnIID(20, dnNormal(0,1))
omega ~ dnUniform(0,1)
# The FMutSel0 model from Yang and Nielsen (2008)
Q1 := fnFMutSel0(fnGTR(er, nuc_pi), F, omega)

# The same -- fMutSel0 = GTR(er,nuc_pi) + X3 + MutSel(F) + dNdS(omega)
Q2 := fndNdS( fnMutSelAA( fnX3( fnGTR(er, nuc_pi)), F), omega))");
	help_strings[string("fnFMutSel0")][string("name")] = string(R"(fnFMutSel0)");
	help_references[string("fnFMutSel0")].push_back(RbHelpReference(R"(Yang, Z. and R. Nielsen. Mutation-Selection Models of Codon Substitution and Their Use to Estimate Selective Strengths on Codon Usage.  Mol. Biol. Evol. (2008) 25(3):568--579)",R"(https://doi.org/10.1093/molbev/msm284 )",R"()"));
	help_arrays[string("fnFMutSel0")][string("see_also")].push_back(string(R"(fnCodonGY94, fnCodonMG94, fnFMutSel0, fnMutSel)"));
	help_strings[string("fnFMutSel0")][string("title")] = string(R"(The FMutSel0 model)");
	help_strings[string("fnFreeBinary")][string("name")] = string(R"(fnFreeBinary)");
	help_strings[string("fnFreeK")][string("name")] = string(R"(fnFreeK)");
	help_strings[string("fnFreeSymmetricRateMatrix")][string("name")] = string(R"(fnFreeSymmetricRateMatrix)");
	help_strings[string("fnGTR")][string("description")] = string(R"(The GTR rate matrix.)");
	help_strings[string("fnGTR")][string("details")] = string(R"(The general time-reversible rate matrix:
  Q(i,j) = exchangeRates(i,j) * baseFrequencies[j]
The exchangeRates are symmetric.)");
	help_strings[string("fnGTR")][string("example")] = string(R"(er ~ dnDirichlet( v(1,1,1,1,1,1) )
pi ~ dnDirichlet( v(1,1,1,1) )
Q := fnGTR(er,pi))");
	help_strings[string("fnGTR")][string("name")] = string(R"(fnGTR)");
	help_references[string("fnGTR")].push_back(RbHelpReference(R"(Tavare, S. Some Probabilistic and Statistical Problems in the Analysis of DNA Sequences.  Lectures on Mathematics in the Life Sciences (1986). 17: 57-86)",R"()",R"(http://www.damtp.cam.ac.uk/user/st321/CV_&_Publications_files/STpapers-pdf/T86.pdf)"));
	help_strings[string("fnGTR")][string("title")] = string(R"(The General Time-Reversible rate matrix)");
	help_strings[string("fnGeographicalDistance")][string("name")] = string(R"(fnGeographicalDistance)");
	help_strings[string("fnHKY")][string("description")] = string(R"(The HKY85 model.)");
	help_strings[string("fnHKY")][string("example")] = string(R"(kappa ~ dnLognormal(0,1)
pi ~ dnDirichlet( v(1,1,1,1) )
Q := fnHKY(kappa,pi))");
	help_strings[string("fnHKY")][string("name")] = string(R"(fnHKY)");
	help_references[string("fnHKY")].push_back(RbHelpReference(R"(Hasegawa, M. et al. Dating of the human-ape splitting by a molecular clock of mitochondrial DNA. Journal of molecular evolution (1985) 22 (2): 160-174.)",R"(https://doi.org/10.1007/BF02101694)",R"(https://link.springer.com/article/10.1007%2FBF02101694 )"));
	help_strings[string("fnHKY")][string("title")] = string(R"(The Hasegawa-Kishino-Yano (1985) nucleotide rate matrix)");
	help_strings[string("fnHiddenStateRateMatrix")][string("name")] = string(R"(fnHiddenStateRateMatrix)");
	help_strings[string("fnHostSwitchRateModifier")][string("name")] = string(R"(fnHostSwitchRateModifier)");
	help_strings[string("fnInfiniteSites")][string("name")] = string(R"(fnInfiniteSites)");
	help_strings[string("fnJC")][string("name")] = string(R"(fnJC)");
	help_strings[string("fnJones")][string("name")] = string(R"(fnJones)");
	help_strings[string("fnK80")][string("name")] = string(R"(fnK80)");
	help_strings[string("fnK81")][string("name")] = string(R"(fnK81)");
	help_strings[string("fnLG")][string("name")] = string(R"(fnLG)");
	help_strings[string("fnLnProbability")][string("name")] = string(R"(fnLnProbability)");
	help_strings[string("fnMixtureCladoProbs")][string("name")] = string(R"(fnMixtureCladoProbs)");
	help_strings[string("fnMtMam")][string("name")] = string(R"(fnMtMam)");
	help_strings[string("fnMtRev")][string("name")] = string(R"(fnMtRev)");
	help_strings[string("fnMutSel")][string("description")] = string(R"(Constructs a rate matrix from scaled selection coefficients w[i] and
mutation rate matrix mu(i,j).

fnMutSel takes 61 scaled selection coefficients, one for each codon.
This differs from fnMutSelAA, which takes 20 scaled selection coefficients,
one for each amino acid.

A substitution from allele i -> j can be decomposed into
 (1) all individuals initially have state i
 (2) a single individual mutates from i -> j, at rate mu(i,j)
 (3) the allele j goes to fixation

Then the substitution rate Q is then given by
  Q(i,j) = mu(i,j) * Pr(j goes to fixation | i was fixed previously).

The probability of fixation is determined by scaled selection coefficients:
  F[i] = 2*N*s[i]
and the initial frequency 1/N of allele j.)");
	help_strings[string("fnMutSel")][string("example")] = string(R"(er ~ dnDirichlet( v(1,1,1,1,1,1) )
nuc_pi ~ dnDirichlet( rep(2.0, 4) )
F ~ dnIID(61, dnNormal(0,1))
Q := fnMutSel(fnX3(fnGTR(er, nuc_pi) ), F)       # GTR + X3 + MutSel

# A mutation-selection balance model on RNA, with GTR mutation.
F2 ~ dnIID(16, dnNormal(0,1))
Q2 := fnMutSel(fnX2(fnGTR(er,nuc_pi) ), F2)      # GTR + X2 + MutSel)");
	help_strings[string("fnMutSel")][string("name")] = string(R"(fnMutSel)");
	help_references[string("fnMutSel")].push_back(RbHelpReference(R"(Yang, Z. and R. Nielsen. Mutation-Selection Models of Codon Substitution and Their Use to Estimate Selective Strengths on Codon Usage.  Mol. Biol. Evol. (2008) 25(3):568--579)",R"(https://doi.org/10.1093/molbev/msm284 )",R"()"));
	help_arrays[string("fnMutSel")][string("see_also")].push_back(string(R"(fnCodonGY94, fnCodonMG94, fnMutSelAA, fnFMutSel, fndNdS)"));
	help_strings[string("fnMutSel")][string("title")] = string(R"(Add mutation-selection balance to a rate matrix.)");
	help_strings[string("fnMutSelAA")][string("description")] = string(R"(Constructs a rate matrix from scaled selection coefficients w[i] and
mutation rate matrix mu(i,j).

fnMutSelAA takes 20 scaled selection coefficients, one for each amino acid.
This differs from fnMutSel, which takes 61 scaled selection coefficients,
one for each codon.  fnMutSelAA assumes that codons for the same amino acid
have the same fitness.

A substitution from allele i -> j can be decomposed into
 (1) all individuals initially have state i
 (2) a single individual mutates from i -> j, at rate mu(i,j)
 (3) the allele j goes to fixation

Then the substitution rate Q is then given by
  Q(i,j) = mu(i,j) * Pr(j goes to fixation | i was fixed previously).

The probability of fixation is determined by scaled selection coefficients:
  F[i] = 2*N*s[i]
and the initial frequency 1/N of allele j.)");
	help_strings[string("fnMutSelAA")][string("example")] = string(R"(er ~ dnDirichlet( v(1,1,1,1,1,1) )
nuc_pi ~ dnDirichlet( rep(2.0, 4) )
F ~ dnIID(20, dnNormal(0,1))
Q := fnMutSelAA(fnX3(fnGTR(er, nuc_pi)), F))");
	help_strings[string("fnMutSelAA")][string("name")] = string(R"(fnMutSelAA)");
	help_references[string("fnMutSelAA")].push_back(RbHelpReference(R"(Yang, Z. and R. Nielsen. Mutation-Selection Models of Codon Substitution and Their Use to Estimate Selective Strengths on Codon Usage.  Mol. Biol. Evol. (2008) 25(3):568--579)",R"(https://doi.org/10.1093/molbev/msm284 )",R"()"));
	help_arrays[string("fnMutSelAA")][string("see_also")].push_back(string(R"(fnCodonGY94, fnCodonMG94, fnX3, fndNdS, fnMutSel)"));
	help_strings[string("fnMutSelAA")][string("title")] = string(R"(Add mutation-selection balance to a rate matrix -- fitnesses on amino acids)");
	help_strings[string("fnNormalizedQuantile")][string("name")] = string(R"(fnNormalizedQuantile)");
	help_strings[string("fnNumUniqueInVector")][string("name")] = string(R"(fnNumUniqueInVector)");
	help_strings[string("fnOrderedRateMatrix")][string("name")] = string(R"(fnOrderedRateMatrix)");
	help_strings[string("fnPD")][string("name")] = string(R"(fnPD)");
	help_strings[string("fnPartialToCorr")][string("name")] = string(R"(fnPartialToCorr)");
	help_strings[string("fnPattersonsD")][string("name")] = string(R"(fnPattersonsD)");
	help_strings[string("fnPhylogeneticIndependentContrasts")][string("name")] = string(R"(fnPhylogeneticIndependentContrasts)");
	help_strings[string("fnPhylogeneticIndependentContrastsMultiSample")][string("name")] = string(R"(fnPhylogeneticIndependentContrastsMultiSample)");
	help_strings[string("fnPoMo")][string("name")] = string(R"(fnPoMo)");
	help_strings[string("fnPruneTree")][string("name")] = string(R"(fnPruneTree)");
	help_strings[string("fnRangeEvolutionRateModifier")][string("name")] = string(R"(fnRangeEvolutionRateModifier)");
	help_strings[string("fnRateGeneratorSequence")][string("name")] = string(R"(fnRateGeneratorSequence)");
	help_strings[string("fnReversiblePoMo")][string("name")] = string(R"(fnReversiblePoMo)");
	help_strings[string("fnRtRev")][string("name")] = string(R"(fnRtRev)");
	help_strings[string("fnSampledCladogenesisRootFrequencies")][string("name")] = string(R"(fnSampledCladogenesisRootFrequencies)");
	help_strings[string("fnSegregatingSites")][string("name")] = string(R"(fnSegregatingSites)");
	help_strings[string("fnShortestDistance")][string("name")] = string(R"(fnShortestDistance)");
	help_strings[string("fnSiteRateModifier")][string("name")] = string(R"(fnSiteRateModifier)");
	help_strings[string("fnStateCountRateModifier")][string("name")] = string(R"(fnStateCountRateModifier)");
	help_strings[string("fnStirling")][string("name")] = string(R"(fnStirling)");
	help_strings[string("fnStitchTree")][string("name")] = string(R"(fnStitchTree)");
	help_strings[string("fnT92")][string("name")] = string(R"(fnT92)");
	help_strings[string("fnTIM")][string("name")] = string(R"(fnTIM)");
	help_strings[string("fnTVM")][string("name")] = string(R"(fnTVM)");
	help_strings[string("fnTajimasD")][string("name")] = string(R"(fnTajimasD)");
	help_strings[string("fnTajimasPi")][string("name")] = string(R"(fnTajimasPi)");
	help_strings[string("fnTrN")][string("description")] = string(R"(The Tamura-Nei nucleotide rate matrix.)");
	help_strings[string("fnTrN")][string("example")] = string(R"(kappaAG ~ dnLognormal(0,1)    # The purine transition rate
kappaCT ~ dnLognormal(0,1)    # The pyrimindine transition rate
pi ~ dnDirichlet( v(1,1,1,1) )
Q := fnTrN(kappaAT, kappaCT, ,pi))");
	help_strings[string("fnTrN")][string("name")] = string(R"(fnTrN)");
	help_references[string("fnTrN")].push_back(RbHelpReference(R"(Tamura, K. and M. Nei. Estimation of the number of nucleotide substitutions in the control region of mitochondrial DNA in humans and chimpanzees. Molecular biology and evolution (1993) 10(3):512-526.)",R"(https://doi.org/10.1093/oxfordjournals.molbev.a040023)",R"(https://academic.oup.com/mbe/article/10/3/512/1016366 )"));
	help_arrays[string("fnTrN")][string("see_also")].push_back(string(R"(fnHKY)"));
	help_strings[string("fnTrN")][string("title")] = string(R"(The Tamura-Nei (1993) nucleotide rate matrix)");
	help_strings[string("fnTreeAssembly")][string("name")] = string(R"(fnTreeAssembly)");
	help_strings[string("fnTreePairwiseDistances")][string("name")] = string(R"(fnTreePairwiseDistances)");
	help_strings[string("fnTreePairwiseNodalDistances")][string("name")] = string(R"(fnTreePairwiseNodalDistances)");
	help_strings[string("fnTreeScale")][string("name")] = string(R"(fnTreeScale)");
	help_strings[string("fnUpperTriangle")][string("name")] = string(R"(fnUpperTriangle)");
	help_strings[string("fnVT")][string("name")] = string(R"(fnVT)");
	help_strings[string("fnVarCovar")][string("name")] = string(R"(fnVarCovar)");
	help_strings[string("fnWAG")][string("name")] = string(R"(fnWAG)");
	help_strings[string("fnWattersonsTheta")][string("name")] = string(R"(fnWattersonsTheta)");
	help_strings[string("fnX2")][string("description")] = string(R"(Constructs a double rate matrix on the 16 nucleotide pairs.

Rates of change from nucleotide i -> j at each doublet position are given by the
nucleotide rate matrix.  The rate of 2 simultaneous changes is 0.

The X3 function can be used to constructor rate matrices on doublets in a
modular fashion.)");
	help_strings[string("fnX2")][string("example")] = string(R"(
kappa ~ dnLognormal(0,1)
nuc_pi ~ dnDirichlet( rep(2.0, 4) )
# Mutation rate matrix on RNA stems
Q1 := fnX2( fnHKY(kappa, nuc_pi) )
F ~ dnIID(16, dnNormal(0,1))
# Add selection to the rate matrix
Q2 := fnMutSel(Q1, F))");
	help_strings[string("fnX2")][string("name")] = string(R"(fnX2)");
	help_arrays[string("fnX2")][string("see_also")].push_back(string(R"(fnX3)"));
	help_strings[string("fnX2")][string("title")] = string(R"(Construct a doublet (16x16) rate matrix from a nucleotide rate matrix.)");
	help_strings[string("fnX3")][string("description")] = string(R"(Constructs a rate matrix on the 61 non-stop codons (in the standard genetic code).

Rates of change from nucleotide i -> j at each codon position are given by the
nucleotide rate matrix.  The rate of 2 or 3 simultaneous changes is 0.

The X3 function can be used to construct other rate matrices in a modular fashion.
For example:
  (i)  MG94  = F81 + X3 + dNdS
  (ii) MG94K = HKY85 + X3 + dNdS)");
	help_strings[string("fnX3")][string("example")] = string(R"(
kappa ~ dnLognormal(0,1)
omega ~ dnUniform(0,1)
nuc_pi ~ dnDirichlet( rep(2.0, 4) )
Q1 := fnCodonMG94K( kappa, omega, nuc_pi )
# This is the same.
Q2 := fndNdS(fnX3(fnHKY(kappa, nuc_pi)), omega)          # HKY + X3 + dNdS, or HKY*3 + dNdS

er ~ dnDirichlet( v(1,1,1,1,1,1) )
Q3 := fnX3(fnGTR(er, nuc_pi))      # GTR + X3, or GTR*3)");
	help_strings[string("fnX3")][string("name")] = string(R"(fnX3)");
	help_references[string("fnX3")].push_back(RbHelpReference(R"(Redelings, BD (2021). BAli-Phy version 3: Model-based co-estimation of Alignment and Phylogeny.  Bioinformatics (2021) 37(10):3032–3034.)",R"(https://doi.org/10.1093/bioinformatics/btab129 )",R"()"));
	help_arrays[string("fnX3")][string("see_also")].push_back(string(R"(fnCodonGY94, fnCodonMG94K, fndNdS)"));
	help_strings[string("fnX3")][string("title")] = string(R"(Construct a codon rate matrix from a nucleotide rate matrix.)");
	help_strings[string("fnassembleContinuousMRF")][string("name")] = string(R"(fnassembleContinuousMRF)");
	help_strings[string("fndNdS")][string("description")] = string(R"(Constructs a rate matrix on the 61 non-stop codons (in the standard genetic code).

   Q(i,j) = Q'(i,j) * omega if aa(i) != aa(j)
                    * 1     if aa(i) == aa(j)

where aa(i) gives the amino acid for codon i in the standard genetic code, and
Q'(i,j) is the input rate matrix on codons.

The dNdS function can be used to construct other rate matrices in a modular fashion.
For example:
  (i)  MG94  = F81 + X3 + dNdS
  (ii) MG94K = HKY85 + X3 + dNdS)");
	help_strings[string("fndNdS")][string("example")] = string(R"(
kappa ~ dnLognormal(0,1)
omega ~ dnUniform(0,1)
nuc_pi ~ dnDirichlet( rep(2.0, 4) )
Q1 := fnCodonMG94K( kappa, omega, nuc_pi )
# This is the same.
Q2 := fndNdS(fnX3(fnHKY(kappa, nuc_pi)), omega)        # HKY + X3 + dNdS,
                                                       #   or HKY*3 + dNdS

er ~ dnDirichlet( v(1,1,1,1,1,1) )
Q3 := fndNdS(fnX3(fnGTR(er, nuc_pi)), omega)         # GTR + X3 + dNdS)");
	help_strings[string("fndNdS")][string("name")] = string(R"(fndNdS)");
	help_references[string("fndNdS")].push_back(RbHelpReference(R"(Redelings, BD (2021). BAli-Phy version 3: Model-based co-estimation of Alignment and Phylogeny.  Bioinformatics (2021) 37(10):3032–3034.)",R"(https://doi.org/10.1093/bioinformatics/btab129)",R"()"));
	help_arrays[string("fndNdS")][string("see_also")].push_back(string(R"(fnCodonGY94, fnCodonMG94K, fnX3, fnMutSel)"));
	help_strings[string("fndNdS")][string("title")] = string(R"(Add a dN/dS factor to a codon rate matrix.)");
	help_strings[string("formatDiscreteCharacterData")][string("name")] = string(R"(formatDiscreteCharacterData)");
	help_strings[string("gamma")][string("name")] = string(R"(gamma)");
	help_arrays[string("getOption")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("getOption")][string("description")] = string(R"(Get a global option for RevBayes.)");
	help_strings[string("getOption")][string("details")] = string(R"(Options are used to personalize RevBayes and are stored on the local machine. Currently this is rather experimental.)");
	help_strings[string("getOption")][string("example")] = string(R"(# compute the absolute value of a real number
getOption("linewidth")

# let us set the linewidth to a new value
setOption("linewidth", 200)

# now let's check what the value is
getOption("linewidth"))");
	help_strings[string("getOption")][string("name")] = string(R"(getOption)");
	help_arrays[string("getOption")][string("see_also")].push_back(string(R"(setOption)"));
	help_strings[string("getOption")][string("title")] = string(R"(Get a global RevBayes option)");
	help_arrays[string("getwd")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("getwd")][string("description")] = string(R"(Get the current working directory which RevBayes uses.)");
	help_strings[string("getwd")][string("example")] = string(R"(# get the current working directory
getwd()

# let us set a new working directory
setwd("~/Desktop")

# check the working directory again
getwd())");
	help_strings[string("getwd")][string("name")] = string(R"(getwd)");
	help_arrays[string("getwd")][string("see_also")].push_back(string(R"(setwd)"));
	help_strings[string("getwd")][string("title")] = string(R"(Get and print the working directory)");
	help_strings[string("help")][string("description")] = string(R"(Provides general or specific help.)");
	help_strings[string("help")][string("example")] = string(R"(# get general help
help()
# get specific help
help("dnNormal"))");
	help_strings[string("help")][string("name")] = string(R"(help)");
	help_strings[string("help")][string("title")] = string(R"(Get help with RevBayes)");
	help_arrays[string("ifelse")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("ifelse")][string("description")] = string(R"(If the expression is true, then the function returns the first value, otherwise the second value.)");
	help_strings[string("ifelse")][string("details")] = string(R"(The ifelse function is important when the value of a variable should deterministically change during an analysis depending on other variables. Standard if-else statements are not dynamically re-evaluated.)");
	help_strings[string("ifelse")][string("example")] = string(R"(a <- 1
b := ifelse( a == 1, 10, -10 )
b

a <- 2
b)");
	help_strings[string("ifelse")][string("name")] = string(R"(ifelse)");
	help_strings[string("ifelse")][string("title")] = string(R"(If-else statement as a function)");
	help_arrays[string("license")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("license")][string("description")] = string(R"(Print the copyright license of RevBayes.)");
	help_strings[string("license")][string("example")] = string(R"(license())");
	help_strings[string("license")][string("name")] = string(R"(license)");
	help_strings[string("license")][string("title")] = string(R"(Copyright license of RevBayes)");
	help_strings[string("listFiles")][string("name")] = string(R"(listFiles)");
	help_arrays[string("listOptions")][string("authors")].push_back(string(R"(Will Freyman)"));
	help_strings[string("listOptions")][string("description")] = string(R"(List all global options for RevBayes.)");
	help_strings[string("listOptions")][string("details")] = string(R"(Options are used to personalize RevBayes and are stored on the local machine. Currently this is rather experimental.)");
	help_strings[string("listOptions")][string("name")] = string(R"(listOptions)");
	help_arrays[string("listOptions")][string("see_also")].push_back(string(R"(setOption)"));
	help_arrays[string("listOptions")][string("see_also")].push_back(string(R"(getOption)"));
	help_strings[string("listOptions")][string("title")] = string(R"(List global RevBayes options)");
	help_arrays[string("ln")][string("authors")].push_back(string(R"(Jordan Koch)"));
	help_strings[string("ln")][string("description")] = string(R"(Returns the natural log of a (positive) value.)");
	help_strings[string("ln")][string("example")] = string(R"(# create a stochastic node with an exponential distribution
x ~ dnExponential(1)

# create a determinstic node that takes the natural log of x
y := ln(x)

# print the values for x and y
x # x has the stochastic value of 2.940149
y # y has the determined value of 1.07846)");
	help_strings[string("ln")][string("name")] = string(R"(ln)");
	help_arrays[string("ln")][string("see_also")].push_back(string(R"(log)"));
	help_strings[string("ln")][string("title")] = string(R"(Natural log function)");
	help_strings[string("log")][string("name")] = string(R"(log)");
	help_strings[string("logistic")][string("name")] = string(R"(logistic)");
	help_arrays[string("ls")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("ls")][string("description")] = string(R"(Show the content of the workspace.)");
	help_strings[string("ls")][string("details")] = string(R"(The list functions shows all the variables in the current workspace. You can also see all the functions available if you use ls(all=TRUE))");
	help_strings[string("ls")][string("example")] = string(R"(# now we have an empty workspace
ls()
# next wee add a variable
a <- 1
# and we can see it
ls())");
	help_strings[string("ls")][string("name")] = string(R"(ls)");
	help_arrays[string("ls")][string("see_also")].push_back(string(R"(clear)"));
	help_arrays[string("ls")][string("see_also")].push_back(string(R"(exists)"));
	help_strings[string("ls")][string("title")] = string(R"(List workspace content)");
	help_arrays[string("mapTree")][string("authors")].push_back(string(R"(Will Freyman)"));
	help_arrays[string("mapTree")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("mapTree")][string("description")] = string(R"(Finds the maximum a posteriori (MAP) topology from a trace of trees and summarizes branch lengths.)");
	help_strings[string("mapTree")][string("example")] = string(R"(# Read in tree trace
tree_trace = readTreeTrace("output/my.trees", burnin=0.25)

# Generate the MAP tree
map_tree = mapTree(trace=tree_trace, file="map.tree"))");
	help_strings[string("mapTree")][string("name")] = string(R"(mapTree)");
	help_arrays[string("mapTree")][string("see_also")].push_back(string(R"(consensusTree)"));
	help_arrays[string("mapTree")][string("see_also")].push_back(string(R"(mccTree)"));
	help_arrays[string("mapTree")][string("see_also")].push_back(string(R"(treeTrace)"));
	help_arrays[string("mapTree")][string("see_also")].push_back(string(R"(readTreeTrace)"));
	help_strings[string("matrix")][string("name")] = string(R"(matrix)");
	help_strings[string("max")][string("description")] = string(R"(Finds the maximum of a vector of numbers.)");
	help_strings[string("max")][string("example")] = string(R"(a = v(1,2,3,4,5)
max(a)
# this will print 5)");
	help_strings[string("max")][string("name")] = string(R"(max)");
	help_arrays[string("max")][string("see_also")].push_back(string(R"(`min`)"));
	help_strings[string("max")][string("title")] = string(R"(Maximum of a set of numbers)");
	help_arrays[string("maxdiff")][string("authors")].push_back(string(R"(Will Pett)"));
	help_strings[string("maxdiff")][string("description")] = string(R"(Finds the maximum difference in clade probabilities between two posterior samples.)");
	help_strings[string("maxdiff")][string("example")] = string(R"(# Read in tree trace
tree_trace = readTreeTrace("output/my.trees", burnin=0.25, nruns=2)

# Compute the maxdiff statistic
maxdiff = maxdiff(traces=tree_trace))");
	help_strings[string("maxdiff")][string("name")] = string(R"(maxdiff)");
	help_arrays[string("maxdiff")][string("see_also")].push_back(string(R"(readTreeTrace)"));
	help_arrays[string("maximumTree")][string("authors")].push_back(string(R"(Bastien Boussau)"));
	help_strings[string("maximumTree")][string("description")] = string(R"(Builds the maximum species tree given several ultrametric gene trees.)");
	help_strings[string("maximumTree")][string("details")] = string(R"(The maximum species tree is a consistent estimate of the species tree under the multispecies coalescent model, if the gene trees are correct and the effective population size constant along the species tree.)");
	help_strings[string("maximumTree")][string("example")] = string(R"(# We simulate a species tree and gene trees and reconstruct a species tree using maximum tree:
# Let’s simulate a species tree with 10 taxa, 5 gene trees, 1 alleles per species:
n_species <- 10
n_genes <- 5
n_alleles <- 2
# we simulate an ultrametric species tree:
# Species names:
for (i in 1:n_species) {
        species[i] <- taxon(taxonName="Species_"+i, speciesName="Species_"+i)
}
spTree ~ dnBirthDeath(lambda=0.3, mu=0.2, rootAge=10, rho=1, samplingStrategy="uniform", condition="nTaxa", taxa=species)
print(spTree)
# let's pick a constant effective population size of 50:
popSize <- 50
# let's simulate gene trees now:
# taxa names:
for (g in 1:n_genes) {
  for (i in 1:n_species) {
    for (j in 1:n_alleles) {
        taxons[g][(i-1)*n_alleles+j] <- taxon(taxonName="Species_"+i, speciesName="Species_"+i)
    }
  }
  geneTrees[g] ~ dnMultiSpeciesCoalescent(speciesTree=spTree, Ne=popSize, taxa=taxons[g])
}
# Let's compute the maximum tree:
recTree <- maximumTree(geneTrees)
print(recTree))");
	help_strings[string("maximumTree")][string("name")] = string(R"(maximumTree)");
	help_references[string("maximumTree")].push_back(RbHelpReference(R"(High-resolution species trees without concatenation. Scott V. Edwards, Liang Liu, and Dennis K. Pearl . PNAS April 3, 2007 vol. 104 no. 14 .)",R"()",R"(http://www.pnas.org/content/104/14/5936.full )"));
	help_references[string("maximumTree")].push_back(RbHelpReference(R"('Maximum tree: a consistent estimator of the species tree. Liu L, Yu L, Pearl DK.  Journal of Mathematical Biology, 2010. Jan;60(1):95-106.')",R"(https://doi.org/10.1007/s00285-009-0260-0)",R"(https://link.springer.com/article/10.1007%2Fs00285-009-0260-0 )"));
	help_strings[string("maximumTree")][string("title")] = string(R"(Maximum tree function to build a species tree.)");
	help_arrays[string("mccTree")][string("authors")].push_back(string(R"(Will Pett)"));
	help_arrays[string("mccTree")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("mccTree")][string("description")] = string(R"(Finds the maximum clade credibility (MCC) topology from a trace of trees and summarizes branch lengths.)");
	help_strings[string("mccTree")][string("example")] = string(R"(# Read in tree trace
tree_trace = readTreeTrace("output/my.trees", burnin=0.25)

# Generate the MCC tree
map_tree = mccTree(trace=tree_trace, file="mcc.tree"))");
	help_strings[string("mccTree")][string("name")] = string(R"(mccTree)");
	help_arrays[string("mccTree")][string("see_also")].push_back(string(R"(consensusTree)"));
	help_arrays[string("mccTree")][string("see_also")].push_back(string(R"(mapTree)"));
	help_arrays[string("mccTree")][string("see_also")].push_back(string(R"(treeTrace)"));
	help_arrays[string("mccTree")][string("see_also")].push_back(string(R"(readTreeTrace)"));
	help_arrays[string("mcmc")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("mcmc")][string("description")] = string(R"(The MCMC analysis object keeps a model and the associated moves and monitors. The object is used to run Markov chain Monte Carlo (MCMC) simulation on the model, using the provided moves, to obtain a sample of the posterior probability distribution. During the analysis, the monitors are responsible for sampling model parameters of interest.)");
	help_strings[string("mcmc")][string("details")] = string(R"(The MCMC analysis object produced by a call to this function keeps copies of the model and the associated moves and monitors. The MCMC analysis object is used to run Markov chain Monte Carlo (MCMC) simulation on the model, using the provided moves, to obtain a sample of the posterior probability distribution. During the analysis, the monitors are responsible for sampling model parameters of interest.)");
	help_strings[string("mcmc")][string("example")] = string(R"(# Create a simple model (unclamped)
a ~ exponential(1)
mymodel = model(a)

# Create a move vector and a monitor vector
moves[1] = mvScale(a, lambda=1.0, weight=1.0)
monitors[1] = mnFile(a,"output/out.log")

# Create an mcmc object
mymcmcObject = mcmc( mymodel, monitors, moves)

# Run a short analysis
mymcmcObject.burnin( generations = 400, tuningInterval = 100)
mymcmcObject.run( generations = 400)

# print the summary of the operators (now tuned)
mymcmcObject.operatorSummary())");
	help_strings[string("mcmc")][string("name")] = string(R"(mcmc)");
	help_references[string("mcmc")].push_back(RbHelpReference(R"(Metropolis N, AW Rosenbluth, MN Rosenbluth, AH Teller, E Teller (1953). Equation of state calculations by fast computing machines. Journal of Chemical Physics, 21:1087-1092.)",R"()",R"()"));
	help_references[string("mcmc")].push_back(RbHelpReference(R"(Hastings WK (1970) Monte Carlo sampling methods using Markov chains and their applications. Biometrika, 57:97-109.)",R"()",R"()"));
	help_arrays[string("mcmc")][string("see_also")].push_back(string(R"(mcmcmc)"));
	help_strings[string("mcmc")][string("title")] = string(R"(MCMC analysis object)");
	help_arrays[string("mcmcmc")][string("authors")].push_back(string(R"(Michael Landis)"));
	help_arrays[string("mcmcmc")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("mcmcmc")][string("description")] = string(R"(The Mcmcmc analysis object keeps a model and the associated moves and monitors. The object is used to run Metropolis Couped Markov chain Monte Carlo (Mcmcmc) simulation on the model, using the provided moves, to obtain a sample of the posterior probability distribution. During the analysis, the monitors are responsible for sampling model parameters of interest.)");
	help_strings[string("mcmcmc")][string("details")] = string(R"(The Mcmcmc analysis object produced by a call to this function keeps copies of the model and the associated moves and monitors. The Mcmcmc analysis object is used to run Markov chain Monte Carlo (Mcmcmc) simulation on the model, using the provided moves, to obtain a sample of the posterior probability distribution. During the analysis, the monitors are responsible for sampling model parameters of interest.)");
	help_strings[string("mcmcmc")][string("example")] = string(R"(# Create a simple model (unclamped)
a ~ exponential(1)
mymodel = model(a)

# Create a move vector and a monitor vector
moves[1] = mvScale(a, lambda=1.0, weight=1.0)
monitors[1] = mnFile(a,"output/out.log")

# Create an mcmcmc object
myMcmcmcObject = mcmcmc( mymodel, monitors, moves, nchains=4, deltaHeat=5)

# Run a short analysis
myMcmcmcObject.burnin( generations = 400, tuningInterval = 100)
myMcmcmcObject.run( generations = 400)

# print the summary of the operators (now tuned)
myMcmcmcObject.operatorSummary())");
	help_strings[string("mcmcmc")][string("name")] = string(R"(mcmcmc)");
	help_references[string("mcmcmc")].push_back(RbHelpReference(R"("Geyer,C.J. (1991) Markov chain Monte Carlo maximum likelihood. In Keramidas  (ed.), Computing Science and Statistics: Proceedings of the 23rd Symposium on  the Interface. Interface Foundation, Fairfax Station, pp. 156\u2013163.")",R"()",R"()"));
	help_references[string("mcmcmc")].push_back(RbHelpReference(R"("Gilks,W.R. and Roberts,G.O. (1996) Strategies for improving MCMC. In  Gilks,W.R., Richardson,S. and Spiegelhalter (eds) Markov chain Monte Carlo in  Practice. Chapman&Hall, London, 89\u2013114.")",R"()",R"()"));
	help_references[string("mcmcmc")].push_back(RbHelpReference(R"(Altekar, G.; Dwarkadas, S.; Huelsenbeck, J. P. & Ronquist, F. Parallel metropolis coupled Markov chain Monte Carlo for Bayesian phylogenetic inference Bioinformatics, Oxford Univ Press, 2004, 20, 407-415.)",R"()",R"()"));
	help_arrays[string("mcmcmc")][string("see_also")].push_back(string(R"(mcmc)"));
	help_strings[string("mcmcmc")][string("title")] = string(R"(Metropolis-Coupled MCMC analysis object)");
	help_strings[string("mean")][string("description")] = string(R"(Finds the arithmetic mean of a vector of numbers.)");
	help_strings[string("mean")][string("details")] = string(R"(The numbers of the vector are summed and divided by the vector length.)");
	help_strings[string("mean")][string("example")] = string(R"(g = v(2,3,5,6,7)
mean(g)
# 4.6)");
	help_strings[string("mean")][string("name")] = string(R"(mean)");
	help_strings[string("mean")][string("title")] = string(R"(Mean of a vector of numbers)");
	help_strings[string("median")][string("description")] = string(R"(Finds the median of a sorted vector of numbers.)");
	help_strings[string("median")][string("details")] = string(R"(The vector is sorted when `median` is used finding the
number of the sorted values with an equal amount of numbers that
are greater than or less than that value. If the length of the vector is even, there will be no such value. In that case, the two are averaged automatically.)");
	help_strings[string("median")][string("example")] = string(R"(a = v(5,3,2,6,8)
median(a)
# 5 is the result
b = v(1,1,2,3,5,8)
median(b)
# 2.5 is the result)");
	help_strings[string("median")][string("name")] = string(R"(median)");
	help_arrays[string("median")][string("see_also")].push_back(string(R"(`mean`)"));
	help_strings[string("median")][string("title")] = string(R"(Median of a set of numbers)");
	help_strings[string("min")][string("description")] = string(R"(Finds the minimum of a vector of numbers.)");
	help_strings[string("min")][string("example")] = string(R"(a = v(0,1,1,2,3,5,8,13)
min(a)
# will print 0)");
	help_strings[string("min")][string("name")] = string(R"(min)");
	help_arrays[string("min")][string("see_also")].push_back(string(R"(`max`)"));
	help_strings[string("min")][string("title")] = string(R"(Minimum of a set of numbers)");
	help_strings[string("mnAncestralState")][string("name")] = string(R"(mnAncestralState)");
	help_strings[string("mnCharHistoryNewick")][string("name")] = string(R"(mnCharHistoryNewick)");
	help_strings[string("mnCharHistoryNhx")][string("name")] = string(R"(mnCharHistoryNhx)");
	help_strings[string("mnCharacterHistorySummary")][string("name")] = string(R"(mnCharacterHistorySummary)");
	help_strings[string("mnExtNewick")][string("name")] = string(R"(mnExtNewick)");
	help_strings[string("mnFile")][string("name")] = string(R"(mnFile)");
	help_strings[string("mnHomeologPhase")][string("name")] = string(R"(mnHomeologPhase)");
	help_strings[string("mnJointConditionalAncestralState")][string("name")] = string(R"(mnJointConditionalAncestralState)");
	help_strings[string("mnModel")][string("name")] = string(R"(mnModel)");
	help_strings[string("mnNexus")][string("name")] = string(R"(mnNexus)");
	help_strings[string("mnProbability")][string("name")] = string(R"(mnProbability)");
	help_strings[string("mnScreen")][string("name")] = string(R"(mnScreen)");
	help_strings[string("mnStochasticBranchRate")][string("name")] = string(R"(mnStochasticBranchRate)");
	help_strings[string("mnStochasticBranchStateTimes")][string("name")] = string(R"(mnStochasticBranchStateTimes)");
	help_strings[string("mnStochasticCharacterMap")][string("name")] = string(R"(mnStochasticCharacterMap)");
	help_strings[string("mnStochasticVariable")][string("name")] = string(R"(mnStochasticVariable)");
	help_strings[string("model")][string("name")] = string(R"(model)");
	help_strings[string("module")][string("name")] = string(R"(module)");
	help_strings[string("mrcaIndex")][string("name")] = string(R"(mrcaIndex)");
	help_strings[string("mvAVMVN")][string("description")] = string(R"(The adaptive variance multivariate-normal proposal of Baele et al. 2017, uses MCMC samples to fit covariance matrix to parameters.

After user-defined waiting time, proposes using covariance matrix epsilon * I + (1 - epsilon) * sigmaSquared * empirical_matrix.

Internally transforms variables based on whether variables are (finitely) bounded, strictly positive, or simplexed.

Non-simplex-valued vector random variables are untransformed.

Add random variables to the move directly (e.g. branch_rates[1], not branch_rates).WARNING: Disabling tuning disables both tuning of proposal variance and learning of empirical covariance matrix.)");
	help_strings[string("mvAVMVN")][string("name")] = string(R"(mvAVMVN)");
	help_strings[string("mvAddRemoveTip")][string("name")] = string(R"(mvAddRemoveTip)");
	help_strings[string("mvBetaProbability")][string("name")] = string(R"(mvBetaProbability)");
	help_strings[string("mvBetaSimplex")][string("name")] = string(R"(mvBetaSimplex)");
	help_strings[string("mvBinarySwitch")][string("name")] = string(R"(mvBinarySwitch)");
	help_strings[string("mvBirthDeathEvent")][string("name")] = string(R"(mvBirthDeathEvent)");
	help_strings[string("mvBirthDeathEventContinuous")][string("name")] = string(R"(mvBirthDeathEventContinuous)");
	help_strings[string("mvBirthDeathFromAgeEvent")][string("name")] = string(R"(mvBirthDeathFromAgeEvent)");
	help_strings[string("mvBranchLengthScale")][string("name")] = string(R"(mvBranchLengthScale)");
	help_strings[string("mvBurstEvent")][string("name")] = string(R"(mvBurstEvent)");
	help_strings[string("mvCharacterHistory")][string("name")] = string(R"(mvCharacterHistory)");
	help_strings[string("mvCollapseExpandFossilBranch")][string("name")] = string(R"(mvCollapseExpandFossilBranch)");
	help_strings[string("mvConjugateInverseWishart")][string("name")] = string(R"(mvConjugateInverseWishart)");
	help_strings[string("mvContinuousCharacterDataSlide")][string("name")] = string(R"(mvContinuousCharacterDataSlide)");
	help_strings[string("mvContinuousEventScale")][string("name")] = string(R"(mvContinuousEventScale)");
	help_arrays[string("mvCorrelationMatrixElementSwap")][string("authors")].push_back(string(R"(Michael R. May)"));
	help_strings[string("mvCorrelationMatrixElementSwap")][string("description")] = string(R"(Swaps elements i and j of the correlation matrix (i != j).)");
	help_strings[string("mvCorrelationMatrixElementSwap")][string("example")] = string(R"(
# we draw a correlation matrix from an LKJ distribution
R ~ dnLKJ(eta=1, dim=5)

# we specify an element swap move
moves[1] = mvCorrelationMatrixElementSwap(R))");
	help_strings[string("mvCorrelationMatrixElementSwap")][string("name")] = string(R"(mvCorrelationMatrixElementSwap)");
	help_strings[string("mvCorrelationMatrixElementSwap")][string("title")] = string(R"(Correlation Matrix element swap move.)");
	help_strings[string("mvCorrelationMatrixRandomWalk")][string("name")] = string(R"(mvCorrelationMatrixRandomWalk)");
	help_arrays[string("mvCorrelationMatrixSingleElementBeta")][string("authors")].push_back(string(R"(Michael R. May)"));
	help_strings[string("mvCorrelationMatrixSingleElementBeta")][string("description")] = string(R"(Beta proposal on a random element of a correlation matrix.)");
	help_strings[string("mvCorrelationMatrixSingleElementBeta")][string("details")] = string(R"(This move chooses a single element of the correlation matrix at random, and draws a proposed value from a Beta distribution centered on the current value (and stretched to range from -1 to 1).)");
	help_strings[string("mvCorrelationMatrixSingleElementBeta")][string("example")] = string(R"(
# we draw a correlation matrix from an LKJ distribution
R ~ dnLKJ(eta=1, dim=5)

# we specify a beta move on the correlation matrix
moves[1] = mvCorrelationMatrixSingleElementBeta(R, alpha=10.0))");
	help_strings[string("mvCorrelationMatrixSingleElementBeta")][string("name")] = string(R"(mvCorrelationMatrixSingleElementBeta)");
	help_arrays[string("mvCorrelationMatrixSingleElementBeta")][string("see_also")].push_back(string(R"(mvCorrelationMatrixSpecificElementBeta)"));
	help_arrays[string("mvCorrelationMatrixSingleElementBeta")][string("see_also")].push_back(string(R"(mvCorrelationMatrixRandomWalk)"));
	help_strings[string("mvCorrelationMatrixSingleElementBeta")][string("title")] = string(R"(Correlation Matrix Beta proposal.)");
	help_strings[string("mvCorrelationMatrixSpecificElementBeta")][string("name")] = string(R"(mvCorrelationMatrixSpecificElementBeta)");
	help_strings[string("mvCorrelationMatrixUpdate")][string("name")] = string(R"(mvCorrelationMatrixUpdate)");
	help_strings[string("mvDPPAllocateAuxGibbs")][string("name")] = string(R"(mvDPPAllocateAuxGibbs)");
	help_strings[string("mvDPPGibbsConcentration")][string("name")] = string(R"(mvDPPGibbsConcentration)");
	help_strings[string("mvDPPValueBetaSimplex")][string("name")] = string(R"(mvDPPValueBetaSimplex)");
	help_strings[string("mvDPPValueScaling")][string("name")] = string(R"(mvDPPValueScaling)");
	help_strings[string("mvDirichletSimplex")][string("name")] = string(R"(mvDirichletSimplex)");
	help_strings[string("mvDiscreteEventCategoryRandomWalk")][string("name")] = string(R"(mvDiscreteEventCategoryRandomWalk)");
	help_strings[string("mvElementSwapSimplex")][string("name")] = string(R"(mvElementSwapSimplex)");
	help_strings[string("mvEllipticalSliceSamplingSimple")][string("name")] = string(R"(mvEllipticalSliceSamplingSimple)");
	help_arrays[string("mvEmpiricalTree")][string("authors")].push_back(string(R"(Will Freyman)"));
	help_arrays[string("mvEmpiricalTree")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_arrays[string("mvEmpiricalTree")][string("authors")].push_back(string(R"(Will Pett)"));
	help_arrays[string("mvEmpiricalTree")][string("authors")].push_back(string(R"(Jiansi Gao)"));
	help_strings[string("mvEmpiricalTree")][string("description")] = string(R"(An MCMC move that operates on empirical tree distributions.)");
	help_strings[string("mvEmpiricalTree")][string("example")] = string(R"(# Read in tree trace
tree_trace = readTreeTrace("output/my.trees", burnin=0.25)

# Create a distribution of trees
tree ~ dnEmpiricalTree(tree_trace)

# Add an MCMC move
moves[1] = mvEmpiricalTree(tree))");
	help_strings[string("mvEmpiricalTree")][string("name")] = string(R"(mvEmpiricalTree)");
	help_arrays[string("mvEmpiricalTree")][string("see_also")].push_back(string(R"(mvEmpiricalTree)"));
	help_arrays[string("mvEmpiricalTree")][string("see_also")].push_back(string(R"(treeTrace)"));
	help_arrays[string("mvEmpiricalTree")][string("see_also")].push_back(string(R"(readTreeTrace)"));
	help_strings[string("mvEventTimeBeta")][string("name")] = string(R"(mvEventTimeBeta)");
	help_strings[string("mvEventTimeSlide")][string("name")] = string(R"(mvEventTimeSlide)");
	help_strings[string("mvFNPR")][string("name")] = string(R"(mvFNPR)");
	help_strings[string("mvGMRFHyperpriorGibbs")][string("name")] = string(R"(mvGMRFHyperpriorGibbs)");
	help_strings[string("mvGMRFUnevenGridHyperpriorGibbs")][string("name")] = string(R"(mvGMRFUnevenGridHyperpriorGibbs)");
	help_strings[string("mvGPR")][string("name")] = string(R"(mvGPR)");
	help_arrays[string("mvGammaScale")][string("authors")].push_back(string(R"(Jeremy M. Brown)"));
	help_strings[string("mvGammaScale")][string("description")] = string(R"(A move to scale a single continuous value by multiplying by a value drawn from a Gamma(lambda,1) distribution. Lambda is the tuning parameter that controls the size of the proposals.)");
	help_strings[string("mvGammaScale")][string("example")] = string(R"(# Here is a simple example for conducting MCMC on the mean and sd of a Normal distribution.

# Uniform(0,1) priors on the mean and sd
mean ~ dnUnif(0,1)
sd ~ dnUnif(0,1)

# Dummy data (will not actually be analyzed)
data <- v(0.4,0.5,0.6)

# Clamping data
for (i in 1:data.size()){ outcomes[i] ~ dnNorm(mean,sd); outcomes[i].clamp(data[i]) }

# Initializing move and monitor counters
mvi = 1
mni = 1

# Adding Gamma scale moves for the mean and sd (THIS MOVE IS HERE)
moves[mvi++] = mvGammaScale(mean)
moves[mvi++] = mvGammaScale(sd)

# Instantiating the model
mymodel = model(outcomes)

# Adding screen monitor for the mean
monitors[mni++] = mnScreen(mean, printgen=1000)

# Creating MCMC object
mymcmc = mcmc(mymodel, moves, monitors)

# Running MCMC under the prior
mymcmc.run(30000,underPrior=TRUE);)");
	help_strings[string("mvGammaScale")][string("name")] = string(R"(mvGammaScale)");
	help_arrays[string("mvGammaScale")][string("see_also")].push_back(string(R"(mvScale)"));
	help_strings[string("mvGibbsDrawCharacterHistory")][string("name")] = string(R"(mvGibbsDrawCharacterHistory)");
	help_strings[string("mvGibbsMixtureAllocation")][string("name")] = string(R"(mvGibbsMixtureAllocation)");
	help_strings[string("mvGraphFlipClique")][string("name")] = string(R"(mvGraphFlipClique)");
	help_strings[string("mvGraphFlipEdge")][string("name")] = string(R"(mvGraphFlipEdge)");
	help_strings[string("mvGraphShiftEdge")][string("name")] = string(R"(mvGraphShiftEdge)");
	help_strings[string("mvHSRFHyperpriorsGibbs")][string("name")] = string(R"(mvHSRFHyperpriorsGibbs)");
	help_strings[string("mvHSRFUnevenGridHyperpriorsGibbs")][string("name")] = string(R"(mvHSRFUnevenGridHyperpriorsGibbs)");
	help_strings[string("mvHomeologPhase")][string("name")] = string(R"(mvHomeologPhase)");
	help_strings[string("mvIndependentTopology")][string("name")] = string(R"(mvIndependentTopology)");
	help_arrays[string("mvLayeredScaleProposal")][string("authors")].push_back(string(R"(Bastien Boussau)"));
	help_strings[string("mvLayeredScaleProposal")][string("description")] = string(R"(Makes a subtree scale move on all subtrees below a given age in the tree. Tree topology is not altered.)");
	help_strings[string("mvLayeredScaleProposal")][string("details")] = string(R"(The tree must be ultrametric.

An age is randomly drawn between the root age and the age of the oldest tip. Then all subtrees below this age are scaled up or down depending on a scaler drawn from an exponential distribution.)");
	help_strings[string("mvLayeredScaleProposal")][string("example")] = string(R"(# We are going to save the trees we simulate in the folder simulatedTrees:
dataFolder = "simulatedTrees/"
# Let’s simulate a species tree with 10 taxa, 2 gene trees, 3 alleles per species:
n_species <- 10
n_genes <- 2
n_alleles <- 3
# we simulate an ultrametric species tree:
# Species names:
for (i in 1:n_species) {
        species[i] <- taxon(taxonName="Species_"+i, speciesName="Species_"+i)
}
spTree ~ dnBirthDeath(lambda=0.3, mu=0.2, rootAge=10, rho=1, samplingStrategy="uniform", condition="nTaxa", taxa=species)
print(spTree)
# let's pick a constant effective population size of 50:
popSize <- 50
# let's simulate gene trees now:
# taxa names:
for (g in 1:n_genes) {
  for (i in 1:n_species) {
    for (j in 1:n_alleles) {
        taxons[g][(i-1)*n_alleles+j] <- taxon(taxonName="Species_"+i+"_"+j, speciesName="Species_"+i)
    }
  }
  geneTrees[g] ~ dnMultiSpeciesCoalescent(speciesTree=spTree, Ne=popSize, taxa=taxons[g])
  print(geneTrees[g])
}
# We can save the species tree and the gene trees:
write(spTree, filename=dataFolder+"speciesTree")
# Saving the gene trees
for (i in 1:(n_genes)) {
  write(geneTrees[i], filename=dataFolder+"geneTree_"+i+".tree")
}
# set my move index
mi = 0
move_species_subtree_scale = mvLayeredScaleProposal( speciesTree=spTree, weight=5 )
for (i in 1:n_genes) {
   move_species_subtree_scale.addGeneTreeVariable( geneTrees[i] )
}
moves[++mi] = move_species_subtree_scale
# We get a handle on our model.
# We can use any node of our model as a handle, here we choose to use the topology.
mymodel = model(spTree)
# Monitors to check the progression of the program
monitors[1] = mnScreen(printgen=10, spTree)
# Here we use a plain MCMC. You could also set nruns=2 for a replicated analysis
# or use mcmcmc with heated chains.
mymcmc = mcmc(mymodel, monitors, moves, nruns=4)
mymcmc.run(generations=1000)
mymcmc.operatorSummary())");
	help_strings[string("mvLayeredScaleProposal")][string("name")] = string(R"(mvLayeredScaleProposal)");
	help_arrays[string("mvLayeredScaleProposal")][string("see_also")].push_back(string(R"(mvSubTreeScale)"));
	help_strings[string("mvLayeredScaleProposal")][string("title")] = string(R"(Rescales all the subtrees below some age.)");
	help_strings[string("mvLevyJump")][string("name")] = string(R"(mvLevyJump)");
	help_strings[string("mvLevyJumpSum")][string("name")] = string(R"(mvLevyJumpSum)");
	help_strings[string("mvMatrixElementScale")][string("name")] = string(R"(mvMatrixElementScale)");
	help_strings[string("mvMatrixElementSlide")][string("name")] = string(R"(mvMatrixElementSlide)");
	help_strings[string("mvMirror")][string("description")] = string(R"(The adaptive mirror (normal) proposal of Thawornwattana et al. 2017, uses MCMC samples to find posterior mean and variance. After user-defined waiting time, proposes moves on opposite side of posterior mean from current location using a normal distribution with the learned posterior standard deviation (scaled by lambda). Before this time, the move uses mu0 as the mean, and lambda as the standard deviation. WARNING: Disabling tuning disables both tuning of proposal variance and learning of empirical mean and variance. To learn the empirical mean and variance without tuning sigma, set adaptOnly=true.)");
	help_strings[string("mvMirror")][string("name")] = string(R"(mvMirror)");
	help_strings[string("mvMirrorMultiplier")][string("description")] = string(R"(The adaptive mirror multiplier (normal) proposal of Thawornwattana et al. 2017, uses MCMC samples to find posterior mean and variance on the log-scale. After user-defined waiting time, proposes moves (on the log-scale) on opposite side of posterior mean from current location using a normal distribution with the learned posterior standard deviation (scaled by lambda). Before this time, the move uses mu0 as the mean, and lambda as the standard deviation. WARNING: Disabling tuning disables both tuning of proposal variance and learning of empirical mean and variance. To learn the empirical mean and variance without tuning sigma, set adaptOnly=true.)");
	help_strings[string("mvMirrorMultiplier")][string("name")] = string(R"(mvMirrorMultiplier)");
	help_strings[string("mvMixtureAllocation")][string("name")] = string(R"(mvMixtureAllocation)");
	help_strings[string("mvMultiValueEventBirthDeath")][string("name")] = string(R"(mvMultiValueEventBirthDeath)");
	help_strings[string("mvMultiValueEventScale")][string("name")] = string(R"(mvMultiValueEventScale)");
	help_strings[string("mvMultiValueEventSlide")][string("name")] = string(R"(mvMultiValueEventSlide)");
	help_strings[string("mvMultipleElementVectorScale")][string("name")] = string(R"(mvMultipleElementVectorScale)");
	help_strings[string("mvNNI")][string("name")] = string(R"(mvNNI)");
	help_strings[string("mvNarrow")][string("name")] = string(R"(mvNarrow)");
	help_strings[string("mvNarrowExchangeRateMatrix")][string("name")] = string(R"(mvNarrowExchangeRateMatrix)");
	help_strings[string("mvNodeRateTimeSlideUniform")][string("name")] = string(R"(mvNodeRateTimeSlideUniform)");
	help_strings[string("mvNodeTimeScale")][string("name")] = string(R"(mvNodeTimeScale)");
	help_strings[string("mvNodeTimeSlideBeta")][string("name")] = string(R"(mvNodeTimeSlideBeta)");
	help_strings[string("mvNodeTimeSlidePathTruncatedNormal")][string("name")] = string(R"(mvNodeTimeSlidePathTruncatedNormal)");
	help_strings[string("mvNodeTimeSlideUniform")][string("name")] = string(R"(mvNodeTimeSlideUniform)");
	help_strings[string("mvNodeTimeSlideUniformAgeConstrained")][string("name")] = string(R"(mvNodeTimeSlideUniformAgeConstrained)");
	help_strings[string("mvRJSwitch")][string("name")] = string(R"(mvRJSwitch)");
	help_strings[string("mvRandomDive")][string("description")] = string(R"(The multiplicative proposal of Dutta 2012, allows for long-distance moves.

Useful for fat-tailed distributions, possibly for bimoodal distributions.

Variables on [0,infinity) are log-transformed for proposals.)");
	help_strings[string("mvRandomDive")][string("name")] = string(R"(mvRandomDive)");
	help_strings[string("mvRandomGeometricWalk")][string("name")] = string(R"(mvRandomGeometricWalk)");
	help_strings[string("mvRandomIntegerWalk")][string("name")] = string(R"(mvRandomIntegerWalk)");
	help_strings[string("mvRateAgeBetaShift")][string("name")] = string(R"(mvRateAgeBetaShift)");
	help_arrays[string("mvResampleFBD")][string("authors")].push_back(string(R"(Walker Pett)"));
	help_strings[string("mvResampleFBD")][string("description")] = string(R"(This move resamples an oldest occurrence age for a random species in a fossilized birth death process described by `dnFBDRP` or `dnFBDRMatrix`)");
	help_strings[string("mvResampleFBD")][string("details")] = string(R"(Under the hood, FBD fossil data is augmented with oldest occurrence ages for each species, which are automatically marginalized during when the model is sampled using MCMC. These ages can also be resampled manually using this move.)");
	help_strings[string("mvResampleFBD")][string("example")] = string(R"(bd ~ dnFBDRP(lambda=lambda, mu=mu, psi=psi, rho=1, taxa=taxa, resample=FALSE)

moves.append( mvResampleFBD(bd, weight=taxa.size()) ))");
	help_strings[string("mvResampleFBD")][string("name")] = string(R"(mvResampleFBD)");
	help_references[string("mvResampleFBD")].push_back(RbHelpReference(R"(The fossilized birth-death model for the analysis of stratigraphic range data under different speciation modes. Stadler, Tanja et al. Journal of theoretical biology, 447:41-55.)",R"()",R"(https://www.sciencedirect.com/science/article/pii/S002251931830119X )"));
	help_arrays[string("mvResampleFBD")][string("see_also")].push_back(string(R"(dnFossilizedBirthDeathRange)"));
	help_arrays[string("mvResampleFBD")][string("see_also")].push_back(string(R"(dnFossilizedBirthDeathRangeMatrix)"));
	help_strings[string("mvRootTimeScaleBactrian")][string("name")] = string(R"(mvRootTimeScaleBactrian)");
	help_strings[string("mvRootTimeSlideUniform")][string("name")] = string(R"(mvRootTimeSlideUniform)");
	help_strings[string("mvSPR")][string("name")] = string(R"(mvSPR)");
	help_strings[string("mvScale")][string("name")] = string(R"(mvScale)");
	help_strings[string("mvScaleBactrian")][string("name")] = string(R"(mvScaleBactrian)");
	help_strings[string("mvScaleBactrianCauchy")][string("name")] = string(R"(mvScaleBactrianCauchy)");
	help_strings[string("mvShrinkExpand")][string("name")] = string(R"(mvShrinkExpand)");
	help_strings[string("mvShrinkExpandScale")][string("name")] = string(R"(mvShrinkExpandScale)");
	help_strings[string("mvSlice")][string("name")] = string(R"(mvSlice)");
	help_strings[string("mvSlide")][string("name")] = string(R"(mvSlide)");
	help_strings[string("mvSlideBactrian")][string("name")] = string(R"(mvSlideBactrian)");
	help_arrays[string("mvSpeciesNarrow")][string("authors")].push_back(string(R"(Sebastian Hoehna, Bastien Boussau)"));
	help_strings[string("mvSpeciesNarrow")][string("description")] = string(R"(Makes a narrow-exchange move both in the species tree and in the gene trees that contain nodes of the relevant populations.)");
	help_strings[string("mvSpeciesNarrow")][string("details")] = string(R"(The species tree must be ultrametric.

All the gene trees that evolved along the species tree according to some form of multispecies coalescent must be added to the move using the addGeneTreeVariable method.

This move jointly performs narrow exchange moves (Nearest-Neighbor Interchanges without branch length alterations) on the species tree and on gene trees, all of which must be ultrametric.)");
	help_strings[string("mvSpeciesNarrow")][string("example")] = string(R"(# We are going to save the trees we simulate in the folder simulatedTrees:
dataFolder = "simulatedTrees/"
# Let’s simulate a species tree with 10 taxa, 2 gene trees, 3 alleles per species:
n_species <- 10
n_genes <- 2
n_alleles <- 3
# we simulate an ultrametric species tree:
# Species names:
for (i in 1:n_species) {
        species[i] <- taxon(taxonName="Species_"+i, speciesName="Species_"+i)
}
spTree ~ dnBirthDeath(lambda=0.3, mu=0.2, rootAge=10, rho=1, samplingStrategy="uniform", condition="nTaxa", taxa=species)
print(spTree)
# let's pick a constant effective population size of 50:
popSize <- 50
# let's simulate gene trees now:
# taxa names:
for (g in 1:n_genes) {
  for (i in 1:n_species) {
    for (j in 1:n_alleles) {
        taxons[g][(i-1)*n_alleles+j] <- taxon(taxonName="Species_"+i+"_"+j, speciesName="Species_"+i)
    }
  }
  geneTrees[g] ~ dnMultiSpeciesCoalescent(speciesTree=spTree, Ne=popSize, taxa=taxons[g])
  print(geneTrees[g])
}
# We can save the species tree and the gene trees:
write(spTree, filename=dataFolder+"speciesTree")
# Saving the gene trees
for (i in 1:(n_genes)) {
  write(geneTrees[i], filename=dataFolder+"geneTree_"+i+".tree")
}
# set my move index
mi = 0
move_species_narrow_exchange = mvSpeciesNarrow( speciesTree=spTree, weight=5 )
for (i in 1:n_genes) {
   move_species_narrow_exchange.addGeneTreeVariable( geneTrees[i] )
}
moves[++mi] = move_species_narrow_exchange
# We get a handle on our model.
# We can use any node of our model as a handle, here we choose to use the topology.
mymodel = model(spTree)
# Monitors to check the progression of the program
monitors[1] = mnScreen(printgen=10, spTree)
# Here we use a plain MCMC. You could also set nruns=2 for a replicated analysis
# or use mcmcmc with heated chains.
mymcmc = mcmc(mymodel, monitors, moves, nruns=4)
mymcmc.run(generations=1000)
mymcmc.operatorSummary())");
	help_strings[string("mvSpeciesNarrow")][string("name")] = string(R"(mvSpeciesNarrow)");
	help_references[string("mvSpeciesNarrow")].push_back(RbHelpReference(R"("Guided Tree Topology Proposals for Bayesian Phylogenetic Inference. Sebastian  H\xF6hna, Alexei J. Drummond. Syst Biol (2012) 61 (1): 1-11.")",R"(https://doi.org/10.1093/sysbio/syr074)",R"(https://academic.oup.com/sysbio/article-lookup/doi/10.1093/sysbio/syr074 )"));
	help_references[string("mvSpeciesNarrow")].push_back(RbHelpReference(R"(Algorithmic improvements to species delimitation and phylogeny estimation under the multispecies coalescent. Graham Jones.  Journal of Mathematical Biology, 2016.)",R"(https://doi.org/10.1007/s00285-016-1034-0)",R"(http://link.springer.com/article/10.1007/s00285-016-1034-0 )"));
	help_arrays[string("mvSpeciesNarrow")][string("see_also")].push_back(string(R"(mvSpeciesSubtreeScale)"));
	help_arrays[string("mvSpeciesNarrow")][string("see_also")].push_back(string(R"(mvSpeciesSubtreeScaleBeta)"));
	help_arrays[string("mvSpeciesNarrow")][string("see_also")].push_back(string(R"(mvSpeciesNodeTimeSlideUniform)"));
	help_arrays[string("mvSpeciesNarrow")][string("see_also")].push_back(string(R"(mvSpeciesTreeScale)"));
	help_strings[string("mvSpeciesNarrow")][string("title")] = string(R"(Narrow-exchange joint move on species tree and gene trees for multispecies coalescent models.)");
	help_arrays[string("mvSpeciesNodeTimeSlideUniform")][string("authors")].push_back(string(R"(Sebastian Hoehna, Bastien Boussau)"));
	help_strings[string("mvSpeciesNodeTimeSlideUniform")][string("description")] = string(R"(Makes a node time slide move both in the species tree and in the gene trees that contain nodes of the relevant populations. Tree topologies are not altered.)");
	help_strings[string("mvSpeciesNodeTimeSlideUniform")][string("details")] = string(R"(The species tree must be ultrametric.

All the gene trees that evolved along the species tree according to some form of multispecies coalescent must be added to the move using the addGeneTreeVariable method.

This move jointly performs node time slides (branch length alterations, keeping the topologies fixed) on the species tree and on gene trees, all of which must be ultrametric.)");
	help_strings[string("mvSpeciesNodeTimeSlideUniform")][string("example")] = string(R"(# We are going to save the trees we simulate in the folder simulatedTrees:
dataFolder = "simulatedTrees"
# Let’s simulate a species tree with 10 taxa, 2 gene trees, 3 alleles per species:
n_species <- 10
n_genes <- 2
n_alleles <- 3
# we simulate an ultrametric species tree:
# Species names:
for (i in 1:n_species) {
        species[i] <- taxon(taxonName="Species_"+i, speciesName="Species_"+i)
}
spTree ~ dnBirthDeath(lambda=0.3, mu=0.2, rootAge=10, rho=1, samplingStrategy="uniform", condition="nTaxa", taxa=species)
print(spTree)
# let's pick a constant effective population size of 50:
popSize <- 50
# let's simulate gene trees now:
# taxa names:
for (g in 1:n_genes) {
  for (i in 1:n_species) {
    for (j in 1:n_alleles) {
        taxons[g][(i-1)*n_alleles+j] <- taxon(taxonName="Species_"+i+"_"+j, speciesName="Species_"+i)
    }
  }
  geneTrees[g] ~ dnMultiSpeciesCoalescent(speciesTree=spTree, Ne=popSize, taxa=taxons[g])
  print(geneTrees[g])
}
# We can save the species tree and the gene trees:
write(spTree, filename=dataFolder+"speciesTree")
# Saving the gene trees
for (i in 1:(n_genes)) {
  write(geneTrees[i], filename=dataFolder+"geneTree_"+i+".tree")
}
# set my move index
mi = 0
move_species_node_time_slide = mvSpeciesNodeTimeSlideUniform( speciesTree=spTree, weight=5 )
for (i in 1:n_genes) {
   move_species_node_time_slide.addGeneTreeVariable( geneTrees[i] )
}
moves[++mi] = move_species_node_time_slide
# We get a handle on our model.
# We can use any node of our model as a handle, here we choose to use the topology.
mymodel = model(spTree)
# Monitors to check the progression of the program
monitors[1] = mnScreen(printgen=10, spTree)
# Here we use a plain MCMC. You could also set nruns=2 for a replicated analysis
# or use mcmcmc with heated chains.
mymcmc = mcmc(mymodel, monitors, moves, nruns=4)
mymcmc.run(generations=1000)
mymcmc.operatorSummary())");
	help_strings[string("mvSpeciesNodeTimeSlideUniform")][string("name")] = string(R"(mvSpeciesNodeTimeSlideUniform)");
	help_references[string("mvSpeciesNodeTimeSlideUniform")].push_back(RbHelpReference(R"("Guided Tree Topology Proposals for Bayesian Phylogenetic Inference. Sebastian  H\xF6hna, Alexei J. Drummond. Syst Biol (2012) 61 (1): 1-11.")",R"(https://doi.org/10.1093/sysbio/syr074)",R"(https://academic.oup.com/sysbio/article-lookup/doi/10.1093/sysbio/syr074 )"));
	help_references[string("mvSpeciesNodeTimeSlideUniform")].push_back(RbHelpReference(R"(Algorithmic improvements to species delimitation and phylogeny estimation under the multispecies coalescent. Graham Jones.  Journal of Mathematical Biology, 2016.)",R"(https://doi.org/10.1007/s00285-016-1034-0)",R"(http://link.springer.com/article/10.1007/s00285-016-1034-0 )"));
	help_arrays[string("mvSpeciesNodeTimeSlideUniform")][string("see_also")].push_back(string(R"(mvSpeciesSubtreeScale)"));
	help_arrays[string("mvSpeciesNodeTimeSlideUniform")][string("see_also")].push_back(string(R"(mvSpeciesSubtreeScaleBeta)"));
	help_arrays[string("mvSpeciesNodeTimeSlideUniform")][string("see_also")].push_back(string(R"(mvSpeciesNarrow)"));
	help_arrays[string("mvSpeciesNodeTimeSlideUniform")][string("see_also")].push_back(string(R"(mvSpeciesTreeScale)"));
	help_strings[string("mvSpeciesNodeTimeSlideUniform")][string("title")] = string(R"(Node time slide joint move on species tree and gene trees for multispecies coalescent models.)");
	help_arrays[string("mvSpeciesSubtreeScale")][string("authors")].push_back(string(R"(Sebastian Hoehna, Bastien Boussau)"));
	help_strings[string("mvSpeciesSubtreeScale")][string("description")] = string(R"(Makes a subtree scale move both in the species tree and in the gene trees that contain nodes of the relevant populations. Tree topologies are not altered.)");
	help_strings[string("mvSpeciesSubtreeScale")][string("details")] = string(R"(The species tree must be ultrametric.

All the gene trees that evolved along the species tree according to some form of multispecies coalescent must be added to the move using the addGeneTreeVariable method.

This move jointly performs a subtree scale move (a whole subtree is scaled up or down, keeping the topology fixed) on the species tree and on gene trees, all of which must be ultrametric.

How this works: we pick a random node which is not the root.
Then, we uniformly pick an age between the parent and the oldest sampled descendant.
The picked subtree is then scaled to this new age.
All gene-trees that are present in the population will be scaled accordingly.)");
	help_strings[string("mvSpeciesSubtreeScale")][string("example")] = string(R"(# We are going to save the trees we simulate in the folder simulatedTrees:
dataFolder = "simulatedTrees/"
# Let’s simulate a species tree with 10 taxa, 2 gene trees, 3 alleles per species:
n_species <- 10
n_genes <- 2
n_alleles <- 3
# we simulate an ultrametric species tree:
# Species names:
for (i in 1:n_species) {
        species[i] <- taxon(taxonName="Species_"+i, speciesName="Species_"+i)
}
spTree ~ dnBirthDeath(lambda=0.3, mu=0.2, rootAge=10, rho=1, samplingStrategy="uniform", condition="nTaxa", taxa=species)
print(spTree)
# let's pick a constant effective population size of 50:
popSize <- 50
# let's simulate gene trees now:
# taxa names:
for (g in 1:n_genes) {
  for (i in 1:n_species) {
    for (j in 1:n_alleles) {
        taxons[g][(i-1)*n_alleles+j] <- taxon(taxonName="Species_"+i+"_"+j, speciesName="Species_"+i)
    }
  }
  geneTrees[g] ~ dnMultiSpeciesCoalescent(speciesTree=spTree, Ne=popSize, taxa=taxons[g])
  print(geneTrees[g])
}
# We can save the species tree and the gene trees:
write(spTree, filename=dataFolder+"speciesTree")
# Saving the gene trees
for (i in 1:(n_genes)) {
  write(geneTrees[i], filename=dataFolder+"geneTree_"+i+".tree")
}
# set my move index
mi = 0
move_species_subtree_scale = mvSpeciesSubtreeScale( speciesTree=spTree, weight=5 )
for (i in 1:n_genes) {
   move_species_subtree_scale.addGeneTreeVariable( geneTrees[i] )
}
moves[++mi] = move_species_subtree_scale
# We get a handle on our model.
# We can use any node of our model as a handle, here we choose to use the topology.
mymodel = model(spTree)
# Monitors to check the progression of the program
monitors[1] = mnScreen(printgen=10, spTree)
# Here we use a plain MCMC. You could also set nruns=2 for a replicated analysis
# or use mcmcmc with heated chains.
mymcmc = mcmc(mymodel, monitors, moves, nruns=4)
mymcmc.run(generations=1000)
mymcmc.operatorSummary())");
	help_strings[string("mvSpeciesSubtreeScale")][string("name")] = string(R"(mvSpeciesSubtreeScale)");
	help_references[string("mvSpeciesSubtreeScale")].push_back(RbHelpReference(R"("Guided Tree Topology Proposals for Bayesian Phylogenetic Inference. Sebastian  H\xF6hna, Alexei J. Drummond. Syst Biol (2012) 61 (1): 1-11.")",R"(https://doi.org/10.1093/sysbio/syr074)",R"(https://academic.oup.com/sysbio/article-lookup/doi/10.1093/sysbio/syr074 )"));
	help_references[string("mvSpeciesSubtreeScale")].push_back(RbHelpReference(R"(Algorithmic improvements to species delimitation and phylogeny estimation under the multispecies coalescent. Graham Jones.  Journal of Mathematical Biology, 2016.)",R"(https://doi.org/10.1007/s00285-016-1034-0)",R"(http://link.springer.com/article/10.1007/s00285-016-1034-0 )"));
	help_arrays[string("mvSpeciesSubtreeScale")][string("see_also")].push_back(string(R"(mvSpeciesNodeTimeSlideUniform)"));
	help_arrays[string("mvSpeciesSubtreeScale")][string("see_also")].push_back(string(R"(mvSpeciesSubtreeScaleBeta)"));
	help_arrays[string("mvSpeciesSubtreeScale")][string("see_also")].push_back(string(R"(mvSpeciesNarrow)"));
	help_arrays[string("mvSpeciesSubtreeScale")][string("see_also")].push_back(string(R"(mvSpeciesTreeScale)"));
	help_strings[string("mvSpeciesSubtreeScale")][string("title")] = string(R"(Subtree scale move on species tree and gene trees for multispecies coalescent models.)");
	help_arrays[string("mvSpeciesSubtreeScaleBeta")][string("authors")].push_back(string(R"(Sebastian Hoehna, Bastien Boussau)"));
	help_strings[string("mvSpeciesSubtreeScaleBeta")][string("description")] = string(R"(Makes a subtree scale move both in the species tree and in the gene trees that contain nodes of the relevant populations. Tree topologies are not altered. Uses a beta distribution to propose a new age value.)");
	help_strings[string("mvSpeciesSubtreeScaleBeta")][string("details")] = string(R"(The species tree must be ultrametric.

All the gene trees that evolved along the species tree according to some form of multispecies coalescent must be added to the move using the addGeneTreeVariable method.

This move jointly performs a subtree scale move (a whole subtree is scaled up or down, keeping the topology fixed) on the species tree and on gene trees, all of which must be ultrametric.

How this works: we pick a random node which is not the root.
Then, we pick a new age between the parent and the oldest sampled descendant according to a beta distribution.
The picked subtree is then scaled to this new age.
All gene-trees that are present in the population will be scaled accordingly.)");
	help_strings[string("mvSpeciesSubtreeScaleBeta")][string("example")] = string(R"(# We are going to save the trees we simulate in the folder simulatedTrees:
dataFolder = "simulatedTrees/"
# Let’s simulate a species tree with 10 taxa, 2 gene trees, 3 alleles per species:
n_species <- 10
n_genes <- 2
n_alleles <- 3
# we simulate an ultrametric species tree:
# Species names:
for (i in 1:n_species) {
        species[i] <- taxon(taxonName="Species_"+i, speciesName="Species_"+i)
}
spTree ~ dnBirthDeath(lambda=0.3, mu=0.2, rootAge=10, rho=1, samplingStrategy="uniform", condition="nTaxa", taxa=species)
print(spTree)
# let's pick a constant effective population size of 50:
popSize <- 50
# let's simulate gene trees now:
# taxa names:
for (g in 1:n_genes) {
  for (i in 1:n_species) {
    for (j in 1:n_alleles) {
        taxons[g][(i-1)*n_alleles+j] <- taxon(taxonName="Species_"+i+"_"+j, speciesName="Species_"+i)
    }
  }
  geneTrees[g] ~ dnMultiSpeciesCoalescent(speciesTree=spTree, Ne=popSize, taxa=taxons[g])
  print(geneTrees[g])
}
# We can save the species tree and the gene trees:
write(spTree, filename=dataFolder+"speciesTree")
# Saving the gene trees
for (i in 1:(n_genes)) {
  write(geneTrees[i], filename=dataFolder+"geneTree_"+i+".tree")
}
# set my move index
mi = 0
move_species_subtree_scale_beta = mvSpeciesSubtreeScaleBeta( speciesTree=spTree, weight=5 )
for (i in 1:n_genes) {
   move_species_subtree_scale_beta.addGeneTreeVariable( geneTrees[i] )
}
moves[++mi] = move_species_subtree_scale_beta
# We get a handle on our model.
# We can use any node of our model as a handle, here we choose to use the topology.
mymodel = model(spTree)
# Monitors to check the progression of the program
monitors[1] = mnScreen(printgen=10, spTree)
# Here we use a plain MCMC. You could also set nruns=2 for a replicated analysis
# or use mcmcmc with heated chains.
mymcmc = mcmc(mymodel, monitors, moves, nruns=4)
mymcmc.run(generations=1000)
mymcmc.operatorSummary())");
	help_strings[string("mvSpeciesSubtreeScaleBeta")][string("name")] = string(R"(mvSpeciesSubtreeScaleBeta)");
	help_references[string("mvSpeciesSubtreeScaleBeta")].push_back(RbHelpReference(R"("Guided Tree Topology Proposals for Bayesian Phylogenetic Inference. Sebastian  H\xF6hna, Alexei J. Drummond. Syst Biol (2012) 61 (1): 1-11.")",R"(https://doi.org/10.1093/sysbio/syr074)",R"(https://academic.oup.com/sysbio/article-lookup/doi/10.1093/sysbio/syr074 )"));
	help_references[string("mvSpeciesSubtreeScaleBeta")].push_back(RbHelpReference(R"(Algorithmic improvements to species delimitation and phylogeny estimation under the multispecies coalescent. Graham Jones.  Journal of Mathematical Biology, 2016.)",R"(https://doi.org/10.1007/s00285-016-1034-0)",R"(http://link.springer.com/article/10.1007/s00285-016-1034-0 )"));
	help_arrays[string("mvSpeciesSubtreeScaleBeta")][string("see_also")].push_back(string(R"(mvSpeciesNodeTimeSlideUniform)"));
	help_arrays[string("mvSpeciesSubtreeScaleBeta")][string("see_also")].push_back(string(R"(mvSpeciesSubtreeScale)"));
	help_arrays[string("mvSpeciesSubtreeScaleBeta")][string("see_also")].push_back(string(R"(mvSpeciesNarrow)"));
	help_arrays[string("mvSpeciesSubtreeScaleBeta")][string("see_also")].push_back(string(R"(mvSpeciesTreeScale)"));
	help_strings[string("mvSpeciesSubtreeScaleBeta")][string("title")] = string(R"(Subtree scale move on species tree and gene trees for multispecies coalescent models.)");
	help_arrays[string("mvSpeciesTreeScale")][string("authors")].push_back(string(R"(Sebastian Hoehna, Bastien Boussau)"));
	help_strings[string("mvSpeciesTreeScale")][string("description")] = string(R"(Makes a tree scale move both in the species tree and in the gene trees. Tree topologies are not altered.)");
	help_strings[string("mvSpeciesTreeScale")][string("details")] = string(R"(The species tree must be ultrametric.

All the gene trees that evolved along the species tree according to some form of multispecies coalescent must be added to the move using the addGeneTreeVariable method.

This move jointly performs a tree scale move (the entire tree is scaled up or down, keeping the topology fixed) on the species tree and on gene trees, all of which must be ultrametric.)");
	help_strings[string("mvSpeciesTreeScale")][string("example")] = string(R"(# We are going to save the trees we simulate in the folder simulatedTrees:
dataFolder = "simulatedTrees/"
# Let’s simulate a species tree with 10 taxa, 2 gene trees, 3 alleles per species:
n_species <- 10
n_genes <- 2
n_alleles <- 3
# we simulate an ultrametric species tree:
# Species names:
for (i in 1:n_species) {
        species[i] <- taxon(taxonName="Species_"+i, speciesName="Species_"+i)
}
root ~  dnNormal(mean=75,sd=2.5,min=0.0, max=Inf)
spTree ~ dnBirthDeath(lambda=0.3, mu=0.2, rootAge=root, rho=1, samplingStrategy="uniform", condition="nTaxa", taxa=species)
print(spTree)
# let's pick a constant effective population size of 50:
popSize <- 50
# let's simulate gene trees now:
# taxa names:
for (g in 1:n_genes) {
  for (i in 1:n_species) {
    for (j in 1:n_alleles) {
        taxons[g][(i-1)*n_alleles+j] <- taxon(taxonName="Species_"+i+"_"+j, speciesName="Species_"+i)
    }
  }
  geneTrees[g] ~ dnMultiSpeciesCoalescent(speciesTree=spTree, Ne=popSize, taxa=taxons[g])
  print(geneTrees[g])
}
# We can save the species tree and the gene trees:
write(spTree, filename=dataFolder+"speciesTree")
# Saving the gene trees
for (i in 1:(n_genes)) {
  write(geneTrees[i], filename=dataFolder+"geneTree_"+i+".tree")
}
# set my move index
mi = 0
move_species_tree_scale = mvSpeciesTreeScale( speciesTree=spTree, root=root, weight=5 )
for (i in 1:n_genes) {
   move_species_tree_scale.addGeneTreeVariable( geneTrees[i] )
}
moves[++mi] = move_species_tree_scale
# We get a handle on our model.
# We can use any node of our model as a handle, here we choose to use the topology.
mymodel = model(spTree)
# Monitors to check the progression of the program
monitors[1] = mnScreen(printgen=10, spTree)
# Here we use a plain MCMC. You could also set nruns=2 for a replicated analysis
# or use mcmcmc with heated chains.
mymcmc = mcmc(mymodel, monitors, moves, nruns=4)
mymcmc.run(generations=1000)
mymcmc.operatorSummary())");
	help_strings[string("mvSpeciesTreeScale")][string("name")] = string(R"(mvSpeciesTreeScale)");
	help_references[string("mvSpeciesTreeScale")].push_back(RbHelpReference(R"("Guided Tree Topology Proposals for Bayesian Phylogenetic Inference. Sebastian  H\xF6hna, Alexei J. Drummond. Syst Biol (2012) 61 (1): 1-11.")",R"(https://doi.org/10.1093/sysbio/syr074)",R"(https://academic.oup.com/sysbio/article-lookup/doi/10.1093/sysbio/syr074 )"));
	help_references[string("mvSpeciesTreeScale")].push_back(RbHelpReference(R"(Algorithmic improvements to species delimitation and phylogeny estimation under the multispecies coalescent. Graham Jones.  Journal of Mathematical Biology, 2016.)",R"(https://doi.org/10.1007/s00285-016-1034-0)",R"(http://link.springer.com/article/10.1007/s00285-016-1034-0 )"));
	help_arrays[string("mvSpeciesTreeScale")][string("see_also")].push_back(string(R"(mvSpeciesNodeTimeSlideUniform)"));
	help_arrays[string("mvSpeciesTreeScale")][string("see_also")].push_back(string(R"(mvSpeciesSubtreeScaleBeta)"));
	help_arrays[string("mvSpeciesTreeScale")][string("see_also")].push_back(string(R"(mvSpeciesNarrow)"));
	help_arrays[string("mvSpeciesTreeScale")][string("see_also")].push_back(string(R"(mvSpeciesSubtreeScale)"));
	help_strings[string("mvSpeciesTreeScale")][string("title")] = string(R"(Tree scale move on species tree and gene trees for multispecies coalescent models.)");
	help_strings[string("mvSubtreeScale")][string("name")] = string(R"(mvSubtreeScale)");
	help_strings[string("mvSymmetricMatrixElementSlide")][string("name")] = string(R"(mvSymmetricMatrixElementSlide)");
	help_strings[string("mvSynchronizedVectorFixedSingleElementSlide")][string("name")] = string(R"(mvSynchronizedVectorFixedSingleElementSlide)");
	help_strings[string("mvTipTimeSlideUniform")][string("name")] = string(R"(mvTipTimeSlideUniform)");
	help_strings[string("mvTreeScale")][string("name")] = string(R"(mvTreeScale)");
	help_strings[string("mvUPPAllocation")][string("name")] = string(R"(mvUPPAllocation)");
	help_strings[string("mvUpDownScale")][string("name")] = string(R"(mvUpDownScale)");
	help_strings[string("mvUpDownSlide")][string("name")] = string(R"(mvUpDownSlide)");
	help_strings[string("mvUpDownSlideBactrian")][string("name")] = string(R"(mvUpDownSlideBactrian)");
	help_strings[string("mvVectorBinarySwitch")][string("name")] = string(R"(mvVectorBinarySwitch)");
	help_strings[string("mvVectorFixedSingleElementSlide")][string("name")] = string(R"(mvVectorFixedSingleElementSlide)");
	help_strings[string("mvVectorScale")][string("name")] = string(R"(mvVectorScale)");
	help_strings[string("mvVectorSingleElementScale")][string("name")] = string(R"(mvVectorSingleElementScale)");
	help_strings[string("mvVectorSingleElementSlide")][string("name")] = string(R"(mvVectorSingleElementSlide)");
	help_strings[string("mvVectorSlide")][string("name")] = string(R"(mvVectorSlide)");
	help_strings[string("mvVectorSlideRecenter")][string("name")] = string(R"(mvVectorSlideRecenter)");
	help_strings[string("nodeAgeByID")][string("name")] = string(R"(nodeAgeByID)");
	help_strings[string("normalize")][string("name")] = string(R"(normalize)");
	help_strings[string("pathSampler")][string("name")] = string(R"(pathSampler)");
	help_strings[string("pomoRF")][string("name")] = string(R"(pomoRF)");
	help_strings[string("pomoState4Convert")][string("name")] = string(R"(pomoState4Convert)");
	help_strings[string("posteriorPredictiveAnalysis")][string("name")] = string(R"(posteriorPredictiveAnalysis)");
	help_strings[string("posteriorPredictiveProbability")][string("name")] = string(R"(posteriorPredictiveProbability)");
	help_strings[string("posteriorPredictiveSimulation")][string("name")] = string(R"(posteriorPredictiveSimulation)");
	help_strings[string("power")][string("name")] = string(R"(power)");
	help_strings[string("powerPosterior")][string("name")] = string(R"(powerPosterior)");
	help_arrays[string("printSeed")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("printSeed")][string("description")] = string(R"(Print the seed of the random number generator.)");
	help_strings[string("printSeed")][string("example")] = string(R"(printSeed()

# Set the seed to a new value
seed(12345)
# Now print the seed again
printSeed())");
	help_strings[string("printSeed")][string("name")] = string(R"(printSeed)");
	help_arrays[string("printSeed")][string("see_also")].push_back(string(R"(seed)"));
	help_strings[string("printSeed")][string("title")] = string(R"(Print the random number generator seed)");
	help_strings[string("quantile")][string("name")] = string(R"(quantile)");
	help_arrays[string("quit")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("quit")][string("description")] = string(R"(Terminates the currently running instance of RevBayes.)");
	help_strings[string("quit")][string("example")] = string(R"(# if you really want to quit
q()
# you can always start again later ...)");
	help_strings[string("quit")][string("name")] = string(R"(quit)");
	help_strings[string("quit")][string("title")] = string(R"(Quit RevBayes)");
	help_arrays[string("range")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("range")][string("description")] = string(R"(Create a sequence of number in the given range (interval).)");
	help_strings[string("range")][string("details")] = string(R"(This function is a simplified version of the sequence function 'seq'. The range function creates a sequence of integer numbers with a step size of 1.)");
	help_strings[string("range")][string("example")] = string(R"(range(1,20)
range(20,-20)

# this function is actually the same as the ':'
20:-20)");
	help_strings[string("range")][string("name")] = string(R"(range)");
	help_arrays[string("range")][string("see_also")].push_back(string(R"(seq)"));
	help_arrays[string("range")][string("see_also")].push_back(string(R"(rep)"));
	help_strings[string("range")][string("title")] = string(R"(A range of consecutive integer numbers)");
	help_strings[string("readAncestralStateTrace")][string("name")] = string(R"(readAncestralStateTrace)");
	help_strings[string("readAncestralStateTreeTrace")][string("name")] = string(R"(readAncestralStateTreeTrace)");
	help_strings[string("readAtlas")][string("name")] = string(R"(readAtlas)");
	help_strings[string("readBranchLengthTrees")][string("name")] = string(R"(readBranchLengthTrees)");
	help_strings[string("readCharacterData")][string("name")] = string(R"(readCharacterData)");
	help_strings[string("readCharacterDataDelimited")][string("name")] = string(R"(readCharacterDataDelimited)");
	help_strings[string("readContinuousCharacterData")][string("name")] = string(R"(readContinuousCharacterData)");
	help_strings[string("readDataDelimitedFile")][string("name")] = string(R"(readDataDelimitedFile)");
	help_strings[string("readDiscreteCharacterData")][string("name")] = string(R"(readDiscreteCharacterData)");
	help_strings[string("readDistanceMatrix")][string("name")] = string(R"(readDistanceMatrix)");
	help_strings[string("readMatrix")][string("name")] = string(R"(readMatrix)");
	help_strings[string("readPoMoCountFile")][string("name")] = string(R"(readPoMoCountFile)");
	help_strings[string("readRelativeNodeAgeConstraints")][string("name")] = string(R"(readRelativeNodeAgeConstraints)");
	help_strings[string("readRelativeNodeAgeWeightedConstraints")][string("name")] = string(R"(readRelativeNodeAgeWeightedConstraints)");
	help_strings[string("readStochasticVariableTrace")][string("name")] = string(R"(readStochasticVariableTrace)");
	help_strings[string("readTaxonData")][string("name")] = string(R"(readTaxonData)");
	help_strings[string("readTrace")][string("name")] = string(R"(readTrace)");
	help_strings[string("readTreeTrace")][string("name")] = string(R"(readTreeTrace)");
	help_arrays[string("readTrees")][string("authors")].push_back(string(R"(Bastien Boussau)"));
	help_strings[string("readTrees")][string("description")] = string(R"(Reads trees from a file containing trees (Nexus, Phylip or Newick accepted), or from a string containing Newick representations of trees.)");
	help_strings[string("readTrees")][string("details")] = string(R"(Either a file name (with the file argument) or a string (with the text argument) must be provided as argument. If both are provided, trees will be read from both sources.)");
	help_strings[string("readTrees")][string("example")] = string(R"(trees=readTrees(text="(a,(b,c));\n(d:0.1,(e:0.1,f:0.1):0.1);")
print(trees))");
	help_strings[string("readTrees")][string("name")] = string(R"(readTrees)");
	help_arrays[string("readTrees")][string("see_also")].push_back(string(R"(write)"));
	help_arrays[string("readTrees")][string("see_also")].push_back(string(R"(readBranchLengthTrees)"));
	help_arrays[string("readTrees")][string("see_also")].push_back(string(R"(readCharacterData)"));
	help_arrays[string("readTrees")][string("see_also")].push_back(string(R"(readCharacterDataDelimited)"));
	help_arrays[string("readTrees")][string("see_also")].push_back(string(R"(readContinuousCharacterData)"));
	help_arrays[string("readTrees")][string("see_also")].push_back(string(R"(readDiscreteCharacterData)"));
	help_arrays[string("readTrees")][string("see_also")].push_back(string(R"(readDataDelimitedFile)"));
	help_arrays[string("readTrees")][string("see_also")].push_back(string(R"(readCharacterData)"));
	help_strings[string("readTrees")][string("title")] = string(R"(Function to read in trees.)");
	help_strings[string("readVCF")][string("name")] = string(R"(readVCF)");
	help_arrays[string("rep")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("rep")][string("description")] = string(R"('rep' creates a vector of 'n' copies of the value 'x'.)");
	help_strings[string("rep")][string("details")] = string(R"('rep' creates a vector of 'n' elements, each with value 'x', preserving the type of 'x' in the returned vector.)");
	help_strings[string("rep")][string("example")] = string(R"(rep(0.1, 3))");
	help_strings[string("rep")][string("name")] = string(R"(rep)");
	help_arrays[string("rep")][string("see_also")].push_back(string(R"(simplex)"));
	help_arrays[string("rep")][string("see_also")].push_back(string(R"(v)"));
	help_strings[string("rep")][string("title")] = string(R"(Replicate a value)");
	help_strings[string("rootedTripletDist")][string("name")] = string(R"(rootedTripletDist)");
	help_strings[string("round")][string("name")] = string(R"(round)");
	help_strings[string("seed")][string("description")] = string(R"(Sets the random number generator seed given a natural number.)");
	help_strings[string("seed")][string("example")] = string(R"(# pick some definitely random number
seed(80797980)
a <- rUniform(1,0.6,1.2)
a
seed(80797980)
a <- rUniform(1,0.6,1.2)
a # this will be the same as above!)");
	help_strings[string("seed")][string("name")] = string(R"(seed)");
	help_strings[string("seed")][string("title")] = string(R"(Seed set function)");
	help_arrays[string("seq")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("seq")][string("description")] = string(R"(Create a sequence of values separate by a given step-size.)");
	help_strings[string("seq")][string("details")] = string(R"(The 'seq' function create a sequence of values, starting with the initial value and then adding the step-size to it until the value reaches the 'to'-value.)");
	help_strings[string("seq")][string("example")] = string(R"(seq(-0.5, 10.5, 2))");
	help_strings[string("seq")][string("name")] = string(R"(seq)");
	help_arrays[string("seq")][string("see_also")].push_back(string(R"(rep)"));
	help_strings[string("seq")][string("title")] = string(R"(Create a sequence values)");
	help_arrays[string("setOption")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("setOption")][string("description")] = string(R"(Set a global option for RevBayes.)");
	help_strings[string("setOption")][string("details")] = string(R"(Options are used to personalize RevBayes and are stored on the local machine. Currently this is rather experimental.)");
	help_strings[string("setOption")][string("example")] = string(R"(# compute the absolute value of a real number
getOption("linewidth")

# let us set the linewidth to a new value
setOption("linewidth", 200)

# now let's check what the value is
getOption("linewidth"))");
	help_strings[string("setOption")][string("name")] = string(R"(setOption)");
	help_arrays[string("setOption")][string("see_also")].push_back(string(R"(getOption)"));
	help_strings[string("setOption")][string("title")] = string(R"(Set a global RevBayes option)");
	help_arrays[string("setwd")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("setwd")][string("description")] = string(R"(Set the current working directory which RevBayes uses.)");
	help_strings[string("setwd")][string("example")] = string(R"(# get the current working directory
getwd()

# let us set a new working directory
setwd("~/Desktop")

# check the working directory again
getwd())");
	help_strings[string("setwd")][string("name")] = string(R"(setwd)");
	help_arrays[string("setwd")][string("see_also")].push_back(string(R"(getwd)"));
	help_strings[string("setwd")][string("title")] = string(R"(Set and print the working directory)");
	help_strings[string("simBirthDeath")][string("description")] = string(R"(Simulates a tree under a very general birth-death process. Parameters are fed in as a n_cats by n_intervals matrix, such that the ith row is the rate vector for the ith category.)");
	help_strings[string("simBirthDeath")][string("name")] = string(R"(simBirthDeath)");
	help_strings[string("simCompleteTree")][string("name")] = string(R"(simCompleteTree)");
	help_strings[string("simStartingTree")][string("name")] = string(R"(simStartingTree)");
	help_strings[string("simTree")][string("name")] = string(R"(simTree)");
	help_strings[string("sinh")][string("name")] = string(R"(sinh)");
	help_strings[string("sort")][string("description")] = string(R"(Function for sorting the members of a vector in either ascending or descending order.)");
	help_strings[string("sort")][string("details")] = string(R"(The vector to be sorted can be of any numeric type. Ascending or descending is specified via the `ascending` argument)");
	help_strings[string("sort")][string("example")] = string(R"(nums = v(1,3,5,7,2,4,6,8)
sort(nums)
# this will result in 1,2,3,4,5,6,7,8
sort(nums, ascending = FALSE)
# this will result in 8,7,6,5,4,3,2,1)");
	help_strings[string("sort")][string("name")] = string(R"(sort)");
	help_strings[string("sort")][string("title")] = string(R"(Sort function)");
	help_strings[string("source")][string("description")] = string(R"(This function takes a Rev filename as an argument and runs that Rev script.)");
	help_strings[string("source")][string("example")] = string(R"(# set the file name
fn = "rb_tutorial.Rev"
# the source function will run the Rev code in the file fn
source(fn))");
	help_strings[string("source")][string("name")] = string(R"(source)");
	help_strings[string("source")][string("title")] = string(R"(Function for sourcing a Rev file)");
	help_strings[string("sqrt")][string("description")] = string(R"(Takes the square root of some positive number `x`.)");
	help_strings[string("sqrt")][string("example")] = string(R"(# compute the square root of a real number
x <- 3.0
root <- sqrt(x)
if ( abs(root*root - x) > 1.0e-15) {
    print("Problem computing the square root.")
} else {
    print("Correct computation of the square root.")
})");
	help_strings[string("sqrt")][string("name")] = string(R"(sqrt)");
	help_arrays[string("sqrt")][string("see_also")].push_back(string(R"(`power`)"));
	help_strings[string("sqrt")][string("title")] = string(R"(Square root of a number)");
	help_strings[string("srGelmanRubin")][string("name")] = string(R"(srGelmanRubin)");
	help_strings[string("srGeweke")][string("name")] = string(R"(srGeweke)");
	help_strings[string("srMaxIteration")][string("name")] = string(R"(srMaxIteration)");
	help_strings[string("srMaxTime")][string("name")] = string(R"(srMaxTime)");
	help_strings[string("srMinESS")][string("name")] = string(R"(srMinESS)");
	help_strings[string("srStationarity")][string("name")] = string(R"(srStationarity)");
	help_strings[string("stdev")][string("name")] = string(R"(stdev)");
	help_strings[string("steppingStoneSampler")][string("name")] = string(R"(steppingStoneSampler)");
	help_arrays[string("stochasticMatrix")][string("authors")].push_back(string(R"(Michael R. May)"));
	help_strings[string("stochasticMatrix")][string("description")] = string(R"(A stochastic matrix is a matrix (not necessarily square) with rows that sum to 1.)");
	help_strings[string("stochasticMatrix")][string("example")] = string(R"(vec[1] ~ dnDirichlet( [1,1,1,1] )
vec[2] ~ dnDirichlet( [1,1,1,1] )
vec[3] ~ dnDirichlet( [1,1,1,1] )
vec[4] ~ dnDirichlet( [1,1,1,1] )

m := stochasticMatrix(vec))");
	help_strings[string("stochasticMatrix")][string("name")] = string(R"(stochasticMatrix)");
	help_strings[string("stochasticMatrix")][string("title")] = string(R"(Building a stochastic matrix.)");
	help_arrays[string("structure")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("structure")][string("description")] = string(R"(Shows all the information about a given variable.)");
	help_strings[string("structure")][string("example")] = string(R"(# create a variable
a <- 1
b := exp(a)
# now create a deterministic variable which will be a child of 'b'
c := ln(b)
# now create a constant variable which will not be a child of 'b'
d <- ln(b)

str(b))");
	help_strings[string("structure")][string("name")] = string(R"(structure)");
	help_arrays[string("structure")][string("see_also")].push_back(string(R"(type)"));
	help_strings[string("structure")][string("title")] = string(R"(The structure of a variable)");
	help_strings[string("sum")][string("description")] = string(R"(Sums all members of a vector of type `Real`, `RealPos`, `Integer`,
or `Natural`)");
	help_strings[string("sum")][string("example")] = string(R"(a = v(1,2,3,4,5,6,7,8)
sum(a)
# returns 36)");
	help_strings[string("sum")][string("name")] = string(R"(sum)");
	help_strings[string("sum")][string("title")] = string(R"(Sum function)");
	help_strings[string("summarizeCharacterMaps")][string("name")] = string(R"(summarizeCharacterMaps)");
	help_strings[string("symmetricDifference")][string("name")] = string(R"(symmetricDifference)");
	help_arrays[string("system")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("system")][string("description")] = string(R"(Run a system command.)");
	help_strings[string("system")][string("details")] = string(R"(This function will delegate the command to the system. In that way, the function works as an interface to the shell.)");
	help_strings[string("system")][string("example")] = string(R"(# We can execute any command just as if you are using a terminal
system("ls")
system("pwd"))");
	help_strings[string("system")][string("name")] = string(R"(system)");
	help_strings[string("system")][string("title")] = string(R"(Execute a system command.)");
	help_strings[string("tanh")][string("name")] = string(R"(tanh)");
	help_arrays[string("taxon")][string("authors")].push_back(string(R"(Michael Landis)"));
	help_strings[string("taxon")][string("description")] = string(R"(The taxon function creates a Taxon object.)");
	help_strings[string("taxon")][string("details")] = string(R"(Each Taxon object records that taxon's name in addition to other information, such as age (which is non-zero for fossils). Character matrices and trees contain Taxon vectors (Taxon[]) that are used to match leaf nodes to data entries for phylogenetic analyses. For multispecies coalescent analyses, Taxon objects are also used to assign species memberships to individuals.)");
	help_strings[string("taxon")][string("example")] = string(R"(# we can create a Taxon object
taxon_gorilla = taxon("Gorilla_gorilla")
# we can create a dummy vector of Taxon objects for simulation
for (i in 1:10) { taxa[i] = taxon("Taxon"+i) }
phy ~ dnBDP(lambda=1, mu=0, rootAge=1, taxa=taxa)
# retrieve the taxon list for 'phy'
phy.taxa())");
	help_strings[string("taxon")][string("name")] = string(R"(taxon)");
	help_arrays[string("taxon")][string("see_also")].push_back(string(R"(readTaxonData)"));
	help_strings[string("taxon")][string("title")] = string(R"(Taxon object)");
	help_arrays[string("time")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("time")][string("description")] = string(R"(Get the current system time in milliseconds.)");
	help_strings[string("time")][string("example")] = string(R"(time()

# Wait a little bit
sum = 0
for (i in 1:10000) sum += i
# Now print the time again
time())");
	help_strings[string("time")][string("name")] = string(R"(time)");
	help_strings[string("time")][string("title")] = string(R"(Get the time information)");
	help_strings[string("tmrca")][string("description")] = string(R"(Finds the most recent common ancestor (TMRCA) of a clade of taxa on a tree.)");
	help_strings[string("tmrca")][string("example")] = string(R"(# let's make up some taxa
taxa = v("horse", "whale", "unicorn", "narwhal")
# convert these to the taxon datatype
for(i in 1:4) { taxa[i] = taxon(taxa[i]) }
# simulate a tree
tau ~ dnUniformTimeTree(rootAge=1, taxa=taxa)
# we also need a molecular substitution model
molecular_model := fnJC(4)
# together these form a continuous time Markov chain over the tree
full_model ~ dnPhyloCTMC(tree=tau, Q=molecular_model, nSites = 100, type="DNA")
# need to make a clade
horned_animals <- clade(taxa[3], taxa[4])
tmrca(tau, horned_animals))");
	help_strings[string("tmrca")][string("name")] = string(R"(tmrca)");
	help_arrays[string("tmrca")][string("see_also")].push_back(string(R"(`clade`)"));
	help_strings[string("tmrca")][string("title")] = string(R"(Find the time to the most recent common ancestor)");
	help_arrays[string("treeTrace")][string("authors")].push_back(string(R"(Will Freyman)"));
	help_strings[string("treeTrace")][string("description")] = string(R"(Creates a tree trace object from a vector of trees.)");
	help_strings[string("treeTrace")][string("example")] = string(R"(# Read in a vector of trees
trees = readTrees("trees.nex")

# Create a tree trace
tree_trace = treeTrace(trees, burnin=0.25)

# Create a distribution of trees from the tree trace
tree ~ dnEmpiricalTree(tree_trace)

# Add an MCMC move
moves[1] = mvEmpiricalTree(tree))");
	help_strings[string("treeTrace")][string("name")] = string(R"(treeTrace)");
	help_arrays[string("treeTrace")][string("see_also")].push_back(string(R"(mvEmpiricalTree)"));
	help_arrays[string("treeTrace")][string("see_also")].push_back(string(R"(treeTrace)"));
	help_arrays[string("treeTrace")][string("see_also")].push_back(string(R"(readTreeTrace)"));
	help_arrays[string("treeTrace")][string("see_also")].push_back(string(R"(readTrees)"));
	help_strings[string("trunc")][string("name")] = string(R"(trunc)");
	help_arrays[string("type")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("type")][string("description")] = string(R"(The value type of a variable.)");
	help_strings[string("type")][string("example")] = string(R"(a <- 2
type(a)

b <- 2.0
type(b))");
	help_strings[string("type")][string("name")] = string(R"(type)");
	help_arrays[string("type")][string("see_also")].push_back(string(R"(structure)"));
	help_strings[string("type")][string("title")] = string(R"(The value type of a variable)");
	help_arrays[string("v")][string("authors")].push_back(string(R"(Michael Landis)"));
	help_strings[string("v")][string("description")] = string(R"('v' creates a vector of the elements '...')");
	help_strings[string("v")][string("details")] = string(R"('v' creates a vector of the elements '...', which are objects of a common base type. Vector elements may themselves be vectors.)");
	help_strings[string("v")][string("example")] = string(R"(# create a vector, Natural[]
x <- v(1,2,3)
x <- x + 1
x

y <- v(2,4,6)
# create a vector of Natural[][]
z <- v(x,y)
z
z[0])");
	help_strings[string("v")][string("name")] = string(R"(v)");
	help_arrays[string("v")][string("see_also")].push_back(string(R"(simplex)"));
	help_arrays[string("v")][string("see_also")].push_back(string(R"(rep)"));
	help_strings[string("v")][string("title")] = string(R"(Create a vector)");
	help_strings[string("validationAnalysis")][string("name")] = string(R"(validationAnalysis)");
	help_strings[string("var")][string("name")] = string(R"(var)");
	help_strings[string("vectorFlatten")][string("name")] = string(R"(vectorFlatten)");
	help_strings[string("write")][string("name")] = string(R"(write)");
	help_strings[string("writeCharacterDataDelimited")][string("name")] = string(R"(writeCharacterDataDelimited)");
	help_strings[string("writeFasta")][string("description")] = string(R"(This function writes out a FASTA formatted file given 
data of class `AbstractHomologousDiscreteCharacterData`.
Filename is specified using the `fn` argument.)");
	help_strings[string("writeFasta")][string("example")] = string(R"(# let's make up some taxa
taxa = v("horse", "whale", "unicorn", "narwhal")
# convert these to the taxon datatype
for(i in 1:4) { taxa[i] = taxon(taxa[i]) }
# simulate a tree
tau ~ dnUniformTimeTree(rootAge=1, taxa=taxa)
# we also need a molecular substitution model
molecular_model := fnJC(4)
# together these form a continuous time Markov chain over the tree
full_model ~ dnPhyloCTMC(tree=tau, Q=molecular_model, nSites = 100, type="DNA")
# this will print a FASTA file with a simulated molecular matrix
# to the working directory
writeFasta(filename="test.fasta", full_model))");
	help_strings[string("writeFasta")][string("name")] = string(R"(writeFasta)");
	help_references[string("writeFasta")].push_back(RbHelpReference(R"(Pearson, William R., and David J. Lipman. "Improved tools for biological sequence comparison." Proceedings of the National Academy of Sciences 85.8 (1988): 2444-2448. )",R"()",R"()"));
	help_references[string("writeFasta")].push_back(RbHelpReference(R"()",R"()",R"(https://www.pnas.org/content/85/8/2444.short )"));
	help_references[string("writeFasta")].push_back(RbHelpReference(R"()",R"(https://doi.org/10.1073/pnas.85.8.2444 )",R"()"));
	help_arrays[string("writeFasta")][string("see_also")].push_back(string(R"(`writeNexus`, `writeCharacterDataDelimited`)"));
	help_strings[string("writeFasta")][string("title")] = string(R"(FASTA file writing function)");
	help_strings[string("writeNexus")][string("description")] = string(R"(Function for writing a nexus file.)");
	help_strings[string("writeNexus")][string("details")] = string(R"(The first argument is the filename  to write to and this must be a string.
The second argument is a data object that must be some character matrix. 
This data matrix could be a morphological matrix, a molecular matrix, or a tree.)");
	help_strings[string("writeNexus")][string("example")] = string(R"(# let's make up some taxa
taxa = v("horse", "whale", "unicorn", "narwhal")
# simulate a tree
tau ~ dnUniformTimeTree(rootAge=1, taxa=taxa)
# we also need a molecular substitution model
molecular_model := fnJC(4)
# together these form a continuous time Markov chain over the tree
full_model ~ dnPhyloCTMC(tree=tau, Q=molecular_model, nSites = 100, type="DNA")
# this will print a Nexus file with a simulated molecular matrix
# to the working directory
writeNexus(filename="test.nex", full_model))");
	help_strings[string("writeNexus")][string("name")] = string(R"(writeNexus)");
	help_references[string("writeNexus")].push_back(RbHelpReference(R"(David R. Maddison, David L. Swofford, Wayne P. Maddison, Nexus: An Extensible File Format for Systematic Information, Systematic Biology, Volume 46, Issue 4, December 1997, Pages 590–621,)",R"(https://doi.org/10.1093/sysbio/46.4.590)",R"(https://academic.oup.com/sysbio/article/46/4/590/1629695)"));
	help_arrays[string("writeNexus")][string("see_also")].push_back(string(R"(`writeFasta`, `writeCharacterDataDelimited`, `write`)"));
	help_strings[string("writeNexus")][string("title")] = string(R"(Nexus file writer)");
}
