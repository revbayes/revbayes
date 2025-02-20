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
	help_strings[string("Simplex")][string("description")] = string(R"(A simplex is a vector of elements that sum to 1.)");
	help_strings[string("Simplex")][string("example")] = string(R"(```rb
x <- simplex([2, 2, 6])
x # = [ 0.2, 0.2, 0.6]
sum(x) # 1, by definition
```)");
	help_strings[string("Simplex")][string("name")] = string(R"(Simplex)");
	help_arrays[string("Simplex")][string("see_also")].push_back(string(R"()"));
	help_arrays[string("Simplex")][string("see_also")].push_back(string(R"(Moves that operate on Simplexes:)"));
	help_arrays[string("Simplex")][string("see_also")].push_back(string(R"(- mvBetaSimplex)"));
	help_arrays[string("Simplex")][string("see_also")].push_back(string(R"(- mvDirichletSimplex)"));
	help_arrays[string("Simplex")][string("see_also")].push_back(string(R"(- mvElementSwapSimplex)"));
	help_strings[string("Simplex")][string("title")] = string(R"(Simplex)");
	help_arrays[string("SiteMixtureModel")][string("authors")].push_back(string(R"(Ben Redelings)"));
	help_strings[string("SiteMixtureModel")][string("description")] = string(R"(A weighted collection of discrete character evolution models.)");
	help_strings[string("SiteMixtureModel")][string("details")] = string(R"(The SiteMixtureModel datatype is a mixture distribution where each
component is a model of discrete character evolution.  Each character evolves
according to one of the component models.  However, the specific model for each
character is not specified in advance.  Instead, each character has some
probability of choosing each component.  These probabilities are specified by
the mixture weights.)");
	help_strings[string("SiteMixtureModel")][string("example")] = string(R"(M := fnInvASRV(fnGammaASRV(fnJC(4),alpha=1),pInv=0.1)
M.weights()
M.nComponents()
M.rootFrequencies(1)

# It possible to express nested models using pipes.
M := fnJC(4) |> fnGammaASRV(alpha=1) |> fnInvASRV(pInv=0.1))");
	help_strings[string("SiteMixtureModel")][string("name")] = string(R"(SiteMixtureModel)");
	help_strings[string("SiteMixtureModel")][string("title")] = string(R"(SiteMixtureModel)");
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
	help_strings[string("clamp")][string("description")] = string(R"(`x.clamp(data)` fixes the value of the stochastic variable `x` to the observation `data`, and marks the variable as corresponding to an observation.)");
	help_strings[string("clamp")][string("details")] = string(R"(Once clamped, the value of `x` is thus expected to remain constant, unless `x` is subsequently unclamped – either explicitly with `x.unclamp()`, or implicitly with `x.clamp(different_data)`.

`x.setValue()` evaluates probabilities at a specific value of `x` without fixing the value.)");
	help_strings[string("clamp")][string("example")] = string(R"(x ~ dnNormal(1, 1)
y ~ dnNormal(2, 2)

# Set the observed value of x
x.clamp(1)
# Compute the probability of the observation
x.probability()

# Modify the observed value of x
x.clamp(2) # equivalent to x.unclamp(); x.clamp(2)
x.probability()

# Evaluate P(y = 1)
y.setValue(1)
y.probability()

# Select another value of y
y.redraw()
print(y)

# Evaluate P(y) at this new value
y.probability()

# Define a model involving x and y
z := x * y

# Because x is clamped, it is invalid to call x.redraw() or mvSlide(x)
# x will remain constant during MCMC, whereas y will be inferred.
mcmc(model(z), [mnScreen(x, y)], [mvSlide(y)]).run(generations = 5))");
	help_strings[string("clamp")][string("name")] = string(R"(Clamp)");
	help_arrays[string("clamp")][string("see_also")].push_back(string(R"(setValue)"));
	help_arrays[string("clamp")][string("see_also")].push_back(string(R"(unclamp)"));
	help_strings[string("clamp")][string("title")] = string(R"(Clamp a stochastic variable to a fixed/observed value)");
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
	help_strings[string("dnBernoulli")][string("description")] = string(R"(The Bernoulli distribution represents a weighted coin toss.)");
	help_strings[string("dnBernoulli")][string("details")] = string(R"(The Bernoulli distribution takes a parameter p, between 0 and 1, and returns 1 with probability p and 0 with probability (1 - p).)");
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
	help_strings[string("dnBeta")][string("description")] = string(R"(The Beta distribution is a flexible distribution that returns a number between 0 and 1, so it is often used as a distribution for probabilities themselves.)");
	help_strings[string("dnBeta")][string("details")] = string(R"(The Beta distribution takes two parameters, alpha and beta. It is equivalent to the uniform when alpha = beta = 1. 

The probability density function is f(x) = x^(alpha - 1) * (1 - x)^(beta - 1) * Gamma(alpha + beta) / (Gamma(alpha) * Gamma(beta)), where Gamma is the gamma function.)");
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
	help_arrays[string("dnBeta")][string("see_also")].push_back(string(R"(gamma)"));
	help_strings[string("dnBeta")][string("title")] = string(R"(Beta Distribution)");
	help_arrays[string("dnBimodalLognormal")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("dnBimodalLognormal")][string("description")] = string(R"(The Bimodal Lognormal distribution unites two separate lognormal distributions.)");
	help_strings[string("dnBimodalLognormal")][string("details")] = string(R"(The bimodal lognormal distribution takes five parameters: mean1 and mean2 (the means (in log-space) of the two lognormal distributions), sd1 and sd2 (the standard deviations (in log-space) of the two lognormal distributions), and p (between 0 and 1). The value will be distributed according to the first lognormal distribution with probability p, and distributed according to the second lognormal distribution with probability (1 - p).)");
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
	help_strings[string("dnBimodalNormal")][string("description")] = string(R"(The Bimodal Normal distribution unites two separate normal distributions.)");
	help_strings[string("dnBimodalNormal")][string("details")] = string(R"(The bimodal normal distribution takes five parameters: mean1 and mean2 (the means of two normal distributions), sd1 and sd2 (the standard deviations of two normal distributions), and p (between 0 and 1). The value will be distributed according to the first normal distribution with probability p, and distributed according to the second normal distribution with probability (1 - p).)");
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
	help_strings[string("dnBimodalNormal")][string("title")] = string(R"(Bimodal Normal dsitribution)");
	help_arrays[string("dnBinomial")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("dnBinomial")][string("description")] = string(R"(The Binomial probability distribution describes the probability of a number of successes for an experiment with a certain number of trials and probability of success per trial.)");
	help_strings[string("dnBinomial")][string("details")] = string(R"(The binomial distribution takes two parameters, p and size. It defines the number of successes in size trials, where each trial has the same success probability p. 

The probability density function is f(x) =  choose(size,x) * p^(x) * (1-p)^(size-p).)");
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
	help_strings[string("dnBivariatePoisson")][string("example")] = string(R"(th0 ~ dnUniform(0.0, 10.0)
    th1 ~ dnUniform(0.0, 10.0)
    th2 ~ dnUniform(0.0, 10.0)
    x ~ dnBivariatePoisson(th0, th1, th2)
    x.clamp([3, 3, 3])
    moves[1] = mvSlide(th0, delta=0.01, weight=1.0)
    moves[2] = mvSlide(th1, delta=0.01, weight=1.0)
    moves[3] = mvSlide(th2, delta=0.01, weight=1.0)
    monitors[1] = mnScreen(printgen=20000, th0)
    mymodel = model(th1)
    mymcmc = mcmc(mymodel, monitors, moves)
    mymcmc.burnin(generations=20000, tuningInterval=100)
    mymcmc.run(generations=200000))");
	help_strings[string("dnBivariatePoisson")][string("name")] = string(R"(dnBivariatePoisson)");
	help_references[string("dnBivariatePoisson")].push_back(RbHelpReference(R"(Karlis D, Ntzoufras J (2003). Bayesian and Non-Bayesian Analysis of Soccer Data using Bivariate Poisson Regression Models. 16th Panhelenic Conference in Statistics, Kavala, April 2003.)",R"()",R"()"));
	help_arrays[string("dnBivariatePoisson")][string("see_also")].push_back(string(R"(dnPoisson)"));
	help_strings[string("dnBivariatePoisson")][string("title")] = string(R"(Bivariate Poisson Distribution)");
	help_strings[string("dnCBDSP")][string("name")] = string(R"(dnCBDSP)");
	help_strings[string("dnCDBDP")][string("name")] = string(R"(dnCDBDP)");
	help_arrays[string("dnCategorical")][string("authors")].push_back(string(R"(Fredrik Ronquist)"));
	help_strings[string("dnCategorical")][string("description")] = string(R"(The Categorical distribution generalizes the Bernoulli distribution, describing the probability of choosing from a number of outcomes, each with their own probability.)");
	help_strings[string("dnCategorical")][string("details")] = string(R"(The categorical distribution takes a parameter p, a simplex (i.e. vector, the elements of which sum to 1). It returns outcome i with probability p[i].

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
	help_strings[string("dnCauchy")][string("description")] = string(R"(The Cauchy distribution describes the distribution of the ratio of two independent normal variables with mean 0 and variance 1.)");
	help_strings[string("dnCauchy")][string("details")] = string(R"(The Cauchy distribution takes two parameters, location and scale. It is a symmetric distribution, but its tails are broad enough that it has no defined mean or variance. The probability density function is f(x) = 1/(pi * scale) * 1 / (1 + x - (location/scale)^2))");
	help_strings[string("dnCauchy")][string("example")] = string(R"(# we simulate some obversations
x <- rCauchy(n=10,location=0,scale=1)
# let's see what the mean and the variance are.
The mean will not converge with more samples, the Cauchy family has no moments.
mean(x)
var(x)
sd(x))");
	help_strings[string("dnCauchy")][string("name")] = string(R"(dnCauchy)");
	help_arrays[string("dnCauchy")][string("see_also")].push_back(string(R"(dnNormal)"));
	help_arrays[string("dnCauchy")][string("see_also")].push_back(string(R"(dnChisq)"));
	help_strings[string("dnCauchy")][string("title")] = string(R"(Cauchy Distribution)");
	help_arrays[string("dnChisq")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("dnChisq")][string("description")] = string(R"(The chi-square distribution with df degrees of freedom describes the distribution of the sum of the squares of df independent normal variables with mean 0 and variance 1.)");
	help_strings[string("dnChisq")][string("details")] = string(R"(The chi-square distribution takes one parameter, df, the number of degrees of freedom. The probability density function is f(x) = x^(df/2 - 1) * e^(-x/2) / (2^(df/2) * Gamma(df/2)), where Gamma is the gamma function.)");
	help_strings[string("dnChisq")][string("example")] = string(R"(# The most important use of the chi-square distribution
# is arguable the quantile function.
# You can access it the following way:
df <- 10
a := qchisq(0.025, df)
a)");
	help_strings[string("dnChisq")][string("name")] = string(R"(dnChisq)");
	help_strings[string("dnChisq")][string("title")] = string(R"(Chi-Square Distribution)");
	help_arrays[string("dnCoalescent")][string("authors")].push_back(string(R"(Ronja Billenstein)"));
	help_arrays[string("dnCoalescent")][string("authors")].push_back(string(R"(Andrew Magee)"));
	help_arrays[string("dnCoalescent")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("dnCoalescent")][string("description")] = string(R"(The constant population size coalescent process specifies a probability density on genealogies, both node ages and the topology.)");
	help_strings[string("dnCoalescent")][string("details")] = string(R"(The underlying theory of the constant population size Coalescent implemented here is Kingman's Coalescent. The implementation here assumes haploid individuals, so for diploid study systems one needs to multiply the effective population size by 2 and the true effective population size in units of individuals needs to be divided by 2 afterwards.
The Coalescent process is parameterized with theta, which here stands for the effective population size (not 4 * Ne * mu). For detailed examples see https://revbayes.github.io/tutorials/coalescent/)");
	help_strings[string("dnCoalescent")][string("example")] = string(R"(
# specify a prior distribution on the constant population size
pop_size ~ dnUniform(0,1E6)
moves.append( mvScale(pop_size, lambda=0.1, tune=true, weight=2.0) )

# specify the coalescent process.
# note that you need to have a vector of taxa
psi ~ dnCoalescent(theta=pop_size, taxa=taxa)

# for monitoring purposes, you may want the root age
root_height := psi.rootAge()

# continue as usual to either clamp the genealogy or infer the genealogy based on sequence data)");
	help_strings[string("dnCoalescent")][string("name")] = string(R"(dnCoalescent)");
	help_references[string("dnCoalescent")].push_back(RbHelpReference(R"(Comparison of Bayesian Coalescent Skyline Plot Models for Inferring Demographic Histories. Billenstein, Ronja and Höhna, Sebastian (2024) Molecular Biology and Evolution, 41(5):msae073.)",R"(https://doi.org/10.1093/molbev/msae073)",R"(https://academic.oup.com/mbe/article/41/5/msae073/7648822 )"));
	help_arrays[string("dnCoalescent")][string("see_also")].push_back(string(R"(dnCoalescentSkyline)"));
	help_arrays[string("dnCoalescent")][string("see_also")].push_back(string(R"(dnCoalescentDemography)"));
	help_strings[string("dnCoalescent")][string("title")] = string(R"(Constant population size Coalescent process)");
	help_strings[string("dnCoalescentDemography")][string("name")] = string(R"(dnCoalescentDemography)");
	help_arrays[string("dnCoalescentSkyline")][string("authors")].push_back(string(R"(Ronja Billenstein)"));
	help_arrays[string("dnCoalescentSkyline")][string("authors")].push_back(string(R"(Andrew Magee)"));
	help_arrays[string("dnCoalescentSkyline")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("dnCoalescentSkyline")][string("description")] = string(R"(The skyline coalescent process specifies a probability density on genealogies, both node ages and the topology. It is used for both heterochronous samples and homochronous samples.)");
	help_strings[string("dnCoalescentSkyline")][string("details")] = string(R"(The underlying theory of the skyline Coalescent implemented here is Kingman's Coalescent. The implementation here assumes haploid individuals, so for diploid study systems one needs to multiply the effective population size by 2 and the true effective population size in units of individuals needs to be divided by 2 afterwards.
The Coalescent process is parameterized with the following parameters:
theta: a vector of effective population sizes (not 4 * Ne * mu).
times: A vector of times for the intervals, if applicable.
events_per_interval: A vector of number of coalescent events for the intervals, if applicable.
method: The method how intervals are defined, either 'specified' or 'events'
model: The shape of the demographic function within the intervals (constant or linear)
taxa: The taxa used when drawing a random tree.
For detailed examples see https://revbayes.github.io/tutorials/coalescent/)");
	help_strings[string("dnCoalescentSkyline")][string("example")] = string(R"(
NUM_INTERVALS = ceil(n_taxa / 5)
for (i in 1:NUM_INTERVALS) {

    pop_size[i] ~ dnUniform(0,1E6)
    pop_size[i].setValue(100.0)
    moves.append( mvScale(pop_size[i], lambda=0.1, tune=true, weight=2.0) )

}

# next we specify a prior on the number of events per interval
# we use a multinomial prior offset to have at least one event per interval
# first, specify the offset
num_events_pi <- rep(1, NUM_INTERVALS)

# next, specify the prior for the multinomial distribution
num_e_simplex_init <- rep(1, NUM_INTERVALS)
num_e_simplex <- simplex(num_e_simplex_init)

# calculate the number of coalescent events that we distribute over the intervals
n_multi <- n_taxa-1-NUM_INTERVALS

# draw the coalescent events into intervals
number_events_pi ~ dnMultinomial(p=num_e_simplex, size=n_multi)

# compute the actual number of events per interval, so the drawn number plus offset
final_number_events_pi := num_events_pi + number_events_pi

moves.append( mvIidPrior(x=number_events_pi) )



### the time tree is a stochastic node modeled by the constant-rate coalescent process (dnCoalescent)
psi ~ dnCoalescentSkyline(theta=pop_size, events_per_interval=final_number_events_pi, method="events", taxa=taxa)

interval_times := psi.getIntervalAges()

root_height := psi.rootAge()


# continue as usual to either clamp the genealogy or infer the genealogy based on sequence data)");
	help_strings[string("dnCoalescentSkyline")][string("name")] = string(R"(dnCoalescentSkyline)");
	help_references[string("dnCoalescentSkyline")].push_back(RbHelpReference(R"(Comparison of Bayesian Coalescent Skyline Plot Models for Inferring Demographic Histories. Billenstein, Ronja and Höhna, Sebastian (2024) Molecular Biology and Evolution, 41(5):msae073.)",R"(https://doi.org/10.1093/molbev/msae073)",R"(https://academic.oup.com/mbe/article/41/5/msae073/7648822 )"));
	help_arrays[string("dnCoalescentSkyline")][string("see_also")].push_back(string(R"(dnCoalescent)"));
	help_arrays[string("dnCoalescentSkyline")][string("see_also")].push_back(string(R"(dnCoalescentDemography)"));
	help_strings[string("dnCoalescentSkyline")][string("title")] = string(R"(Heterochonous and homochronous skyline Coalescent process)");
	help_strings[string("dnCompleteBirthDeath")][string("name")] = string(R"(dnCompleteBirthDeath)");
	help_strings[string("dnConstrainedNodeAge")][string("name")] = string(R"(dnConstrainedNodeAge)");
	help_strings[string("dnConstrainedNodeOrder")][string("name")] = string(R"(dnConstrainedNodeOrder)");
	help_strings[string("dnConstrainedTopology")][string("name")] = string(R"(dnConstrainedTopology)");
	help_strings[string("dnCppNormal")][string("name")] = string(R"(dnCppNormal)");
	help_strings[string("dnDPP")][string("name")] = string(R"(dnDPP)");
	help_strings[string("dnDecomposedInvWishart")][string("name")] = string(R"(dnDecomposedInvWishart)");
	help_arrays[string("dnDirichlet")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("dnDirichlet")][string("description")] = string(R"(The Dirichlet distribution is a generalization of the Beta distribution for multiple variables.)");
	help_strings[string("dnDirichlet")][string("details")] = string(R"(The Dirichlet distribution takes one parameter, alpha, a vector of numbers representing the concentration of the distribution on each variable. It then returns a simplex (i.e. a vector whose elements sum to 1) representing the relative probability of each variable. Note that when every element of alpha is 1, the distribution is equivalent to a uniform on each element.)");
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
	help_strings[string("dnExponential")][string("description")] = string(R"(The Exponential distribution describes the distribution of the times between events in a Poisson point process.)");
	help_strings[string("dnExponential")][string("details")] = string(R"(The exponential distribution takes one parameter, lambda, describing the rate (i.e. 1/mean). The probability density function is f(x) = lambda * exp(-lambda*x).)");
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
	help_strings[string("dnGamma")][string("description")] = string(R"(The Gamma distribution describes the probability of the sum of exponentially distributed variables.)");
	help_strings[string("dnGamma")][string("details")] = string(R"(The gamma distribution takes two parameters, shape and rate. Similar to how 1/rate gives the mean of the exponential, shape/rate gives the mean of the gamma. It provides a natural prior distribution for parameters that could be considered as sums of exponential variables.

The probability density function is f(x) = rate^shape * x^(shape - 1) * e^(-rate * x) / Gamma(shape), where Gamma is the gamma function. Note that for shape = 1, the gamma distribution reduces to an exponential distribution.)");
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
	help_arrays[string("dnLog")][string("authors")].push_back(string(R"(Ben Redelings)"));
	help_strings[string("dnLog")][string("description")] = string(R"(Log-scales a given distribution.)");
	help_strings[string("dnLog")][string("details")] = string(R"(If X ~ dist then exp(X) ~ dnLog(dist)

This provides a way to construct distributions like dnLognormal and
dnLoguniform directly from the underlying distribution in log-space.
It can therefore express distributions that are not directly implemented.

The distribution `dist` can be either univariate (dnNormal) or
multivariate (dnMultivariateNormal).)");
	help_strings[string("dnLog")][string("example")] = string(R"(x ~ dnLog(dnNormal(0,1))          # Draw from the log-Normal distribution
x ~ dnNormal(0,1) |> dnLog()      # Expressed using pipes.
x ~ dnLognormal(0,1)              # This is equivalent.
y ~ dnNormal(0,1)
x := exp(y)                       # This is also equivalent.

x ~ dnLog(dnGamma(2,3))           # There is no equivalent for this.
x ~ dnIID(10,dnLog(dnGamma(2,3))) # Draw 10 log-Gamma(2,3) random variables.

mu = [1.0, 2.0, 3.0, 4.0]
Sigma ~ dnWishart(df=4, kappa=2, dim=4)
x ~ dnLog(dnMultivariateNormal(mu,Sigma)))");
	help_strings[string("dnLog")][string("name")] = string(R"(dnLog)");
	help_arrays[string("dnLog")][string("see_also")].push_back(string(R"(dnLognormal)"));
	help_arrays[string("dnLog")][string("see_also")].push_back(string(R"(dnLoguniform)"));
	help_strings[string("dnLog")][string("title")] = string(R"(Log-scaled distribution)");
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
	help_strings[string("dnPhyloCTMC")][string("description")] = string(R"(dnPhyloCTMC gives the probability distribution of tip labels on a phylogenetic tree given an phylogenetic continuous-time Markov chain model.)");
	help_strings[string("dnPhyloCTMC")][string("details")] = string(R"(
The likelihood of observed tip labels (specified via a clamped `AbstractHomologousDiscreteCharacterData` object) is computed using Felsenstein's pruning algorithm, with partial likelihoods stored for each branch of the tree. It is automatically outputted in the `Likelihood` column of the `mnFile()` and `mnScreen()` monitors (which can be suppressed with `likelihood = FALSE`).)");
	help_strings[string("dnPhyloCTMC")][string("example")] = string(R"(
```rb
# Read character data from a file
chars <- readDiscreteCharacterData("myData.nex")
taxa = chars.taxa()

# Draw a tree with branch lengths
tree ~ dnUniformTopologyBranchLength( taxa, branchLengthDistribution=dnExp(10.0) )

# Define a rate matrix
q_matrix <- fnJC(4)

# Create stochastic node with the tip distribution given by `tree` and `q_matrix`
x ~ dnPhyloCTMC(tree = tree, Q = q_matrix)

# Clamp observed characters to the node
x.clamp(chars)

# Calculate the probability of the observed characters under the given distribution
x.lnProbability()
```)");
	help_strings[string("dnPhyloCTMC")][string("name")] = string(R"(`dnPhyloCTMC`: Distribution of a phylogenetic continuous-time Markov chain)");
	help_references[string("dnPhyloCTMC")].push_back(RbHelpReference(R"(Felsenstein J., 1973. Maximum Likelihood and Minimum-Steps Methods for Estimating Evolutionary Trees from Data on Discrete Characters. Systematic Biology 22:3, 240--249)",R"()",R"()"));
	help_references[string("dnPhyloCTMC")].push_back(RbHelpReference(R"(Felsenstein, J. (1981). Evolutionary trees from DNA sequences: A maximum likelihood approach. Journal of Molecular Evolution. 17 (6): 368–376.)",R"()",R"()"));
	help_references[string("dnPhyloCTMC")].push_back(RbHelpReference(R"(Höhna, S., Landis, M.J. and Heath, T.A. 2017. Phylogenetic inference using `RevBayes`. Curr. Protoc. Bioinform. 57:6.16.1-6.16.34.)",R"(10.1002/cpbi.22)",R"()"));
	help_arrays[string("dnPhyloCTMC")][string("see_also")].push_back(string(R"(- Tutorial on [graphical models](https://revbayes.github.io/tutorials/intro/graph_models))"));
	help_arrays[string("dnPhyloCTMC")][string("see_also")].push_back(string(R"()"));
	help_arrays[string("dnPhyloCTMC")][string("see_also")].push_back(string(R"(- Tutorial on [specifying a phylogenetic continuous-time Markov chain](https://revbayes.github.io/tutorials/ctmc/) model)"));
	help_strings[string("dnPhyloCTMC")][string("title")] = string(R"(The parameters of a phylogenetic model – a tree topology with branch lengths, a substitution model that describes how observations evolve over the tree, etc. – collectively form a distribution called the _phylogenetic continuous-time Markov chain_.)");
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
	help_arrays[string("dnUniformInteger")][string("authors")].push_back(string(R"(Sebastion Hoehna)"));
	help_strings[string("dnUniformInteger")][string("description")] = string(R"(This function creates a stochastic node drawing a random integer from a uniform distribution.)");
	help_strings[string("dnUniformInteger")][string("details")] = string(R"(This function will randomly draw an integer from a uniform distribution
from a minimum and maximum integer set in the first and second arguments.
dnUniformInteger must be defined as a stochastic node. This function can also
be called using the alias 'dnUnifInt'.)");
	help_strings[string("dnUniformInteger")][string("example")] = string(R"(# Create stochastic node
x ~ dnUniformInteger(1, 10)
# See what x was assigned
x
5)");
	help_strings[string("dnUniformInteger")][string("name")] = string(R"(dnUniformInteger)");
	help_arrays[string("dnUniformInteger")][string("see_also")].push_back(string(R"(dnNormal)"));
	help_arrays[string("dnUniformInteger")][string("see_also")].push_back(string(R"(dnExponential)"));
	help_strings[string("dnUniformInteger")][string("title")] = string(R"(Uniform Integer Distribution)");
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
	help_strings[string("exp")][string("description")] = string(R"(Maps the value of a number x to e^x, where e is the number such that ln(e) = 1.)");
	help_strings[string("exp")][string("example")] = string(R"(# checking that ln(e) = 1
x <- exp(1)
ln_of_x <- ln(x)
if (ln_of_ex != 1) {
        print("Problem when computing an exponential value.")
} else {
        print("Correct computation of an exponential value.")
})");
	help_strings[string("exp")][string("name")] = string(R"(exp)");
	help_arrays[string("exp")][string("see_also")].push_back(string(R"(ln)"));
	help_strings[string("exp")][string("title")] = string(R"(Exponential of a number)");
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
	help_strings[string("fnDiscretizeBeta")][string("description")] = string(R"(Select representative values from `num_cats` discrete subdivisions of a beta distribution.)");
	help_strings[string("fnDiscretizeBeta")][string("details")] = string(R"(
A beta distribution is defined by two shape parameters, alpha and beta.

Where a parameter or prior is defined based on the beta distribution, it may be more tractable to evaluate likelihoods at a fixed number of points from the distribution.  These representative points can be computed using `dnDiscretizeBeta`.

In practice, these values are computed as follows:

Let _n_ be the number of categories.
If `median = TRUE`, the quantile function is performed at the midpoint of each category.  Call this vector _q_.
_q_ is then normalized by dividing against its sum, so its elements sum to one; then multiplied by a factor _n_ * _alpha) / (_alpha_ + _beta_).

The computation to obtain the mean for each category, when `median = FALSE`, is more complex, making use of the incomplete beta function ( Majumder & Bhattacharjee 1973).

A real-world use case is available in Wright et al. (2016), with discussion of the properties of the beta distribution. Corresponding tutorials are available at https://www.palass.org/sites/default/files/media/publications/newsletters/number_106/number_106_0.pdf and https://revbayes.github.io/tutorials/morph_tree/V2.html.)");
	help_strings[string("fnDiscretizeBeta")][string("example")] = string(R"(# Values to represent four quadrants of a symmetric beta distribution
categories := fnDiscretizeBeta(0.2, 0.2, 4)
print(categories))");
	help_strings[string("fnDiscretizeBeta")][string("name")] = string(R"(fnDiscretizeBeta)");
	help_references[string("fnDiscretizeBeta")].push_back(RbHelpReference(R"(Majumder & Bhattacharjee. 1973. Algorithm AS63. Applied Statistics, 22.)",R"(NULL)",R"(NULL )"));
	help_references[string("fnDiscretizeBeta")].push_back(RbHelpReference(R"(WRIGHT, A. M., LLOYD, G. T. and HILLIS, D. H. 2016. Modeling character change heterogeneity in phylogenetic analyses of morphology through the use of priors. _Systematic Biology_, 65, 602–11.)",R"(10.1093/sysbio/syv122)",R"(https://doi.org/10.1093/sysbio/syv122 )"));
	help_arrays[string("fnDiscretizeBeta")][string("see_also")].push_back(string(R"(A translation of `fnDiscretizeBeta` into R is available at https://gist.github.com/ms609/883632d10d4d80ea5391cee9c47071fc.)"));
	help_strings[string("fnDiscretizeBeta")][string("title")] = string(R"(Disctetize a beta distribution)");
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
	help_strings[string("fnF81")][string("description")] = string(R"(DNA evolution model proposed in Felsenstein (1981).)");
	help_strings[string("fnF81")][string("details")] = string(R"(In this model, states are allowed to have different stationary frequencies, and exchangeability rates between states are equal. Its only argument, baseFrequencies, codes for said stationary frequencies. While this is usually used for DNA (and therefore has four states), the function can take any number of states, and therefore be used for many other applications (such as aminoacid or morphological evolution).

The F81 rate matrix elements will be of the form:
    Q[i, j] = c * baseFrequencies[j]

where c is a constant needed to normalize the average rate to 1)");
	help_strings[string("fnF81")][string("example")] = string(R"(# stationary base frequencies
baseFrequencies ~ dnDirichlet(v(1,1,1,1))

# create an F81 rate matrix
Q := fnF81(baseFrequencies))");
	help_strings[string("fnF81")][string("name")] = string(R"(fnF81)");
	help_references[string("fnF81")].push_back(RbHelpReference(R"(Felsenstein J (1981). "Evolutionary trees from DNA sequences: a maximum likelihood approach". Journal of Molecular Evolution. 17:368–76.)",R"(https://doi.org/10.1007/BF01734359)",R"(https://link.springer.com/article/10.1007/BF01734359 )"));
	help_arrays[string("fnF81")][string("see_also")].push_back(string(R"(fnJC)"));
	help_arrays[string("fnF81")][string("see_also")].push_back(string(R"(fnK80)"));
	help_arrays[string("fnF81")][string("see_also")].push_back(string(R"(fnK81)"));
	help_arrays[string("fnF81")][string("see_also")].push_back(string(R"(fnT92)"));
	help_arrays[string("fnF81")][string("see_also")].push_back(string(R"(fnHKY)"));
	help_arrays[string("fnF81")][string("see_also")].push_back(string(R"(fnTrN)"));
	help_arrays[string("fnF81")][string("see_also")].push_back(string(R"(fnGTR)"));
	help_strings[string("fnF81")][string("title")] = string(R"(The Felsenstein (1981) rate matrix)");
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
	help_arrays[string("fnFreeK")][string("authors")].push_back(string(R"(Michael Landis)"));
	help_strings[string("fnFreeK")][string("description")] = string(R"(This function generates and returns a free rates matrix.)");
	help_strings[string("fnFreeK")][string("details")] = string(R"(This function accepts both RealPos[] or RealPos[][] as the first argument to automatically
generate a rate matrix with corresponding substitution rates, returning a rate matrix object.
Users can specify if matrix should be normalized in the second argument using a boolean 
variable (default TRUE). Lastly users can specify what matrix exponential method to
use (default eigen) with a string. Possible options include:
scalingAndSquaring
scalingAndSquaringPade
scalingAndSquaringTaylor
uniformization
eigen)");
	help_strings[string("fnFreeK")][string("example")] = string(R"(# Define vector to pass
x <- [1, 1, 1, 1]
# Use fnFreeK to create rate matrix
fnFreeK(x)
[ [ -1.0000, 1.0000 ] ,
  [ 1.0000, -1.0000 ] ])");
	help_strings[string("fnFreeK")][string("name")] = string(R"(fnFreeK)");
	help_arrays[string("fnFreeK")][string("see_also")].push_back(string(R"(RateMatrix)"));
	help_arrays[string("fnFreeK")][string("see_also")].push_back(string(R"(fnFreeBinary)"));
	help_arrays[string("fnFreeK")][string("see_also")].push_back(string(R"(fnFreeSymmetricRateMatrix)"));
	help_strings[string("fnFreeK")][string("title")] = string(R"(Free K Rate Matrix)");
	help_strings[string("fnFreeSymmetricRateMatrix")][string("name")] = string(R"(fnFreeSymmetricRateMatrix)");
	help_strings[string("fnGTR")][string("description")] = string(R"(DNA evolution model proposed in Tavare (1986).)");
	help_strings[string("fnGTR")][string("details")] = string(R"(In this model, states are allowed to have different stationary frequencies, and exchangeability rates between states are allowed to be different. Its first argument, exchangeRates, codes for the transition rates between states (as in other models, transition rates are assumed to be symmetric). Its second argument, baseFrequencies, codes for the stationary frequencies of these states. Note that for n states, exchangeRates should have length n*(n-1)/2, and baseFrequencies should have length n. While this is usually used for DNA (and therefore has four states), the function can take any number of states, and therefore be used for many other applications (such as aminoacid or morphological evolution).

The general time-reversible rate matrix elements will be of the form:
     Q[i, j] = c * exchangeRates[i, j] * baseFrequencies[j]

where c is a constant needed to normalize the average rate to 1.)");
	help_strings[string("fnGTR")][string("example")] = string(R"(# exchange rates
er ~ dnDirichlet( v(1,1,1,1,1,1) )

# base frequencies
pi ~ dnDirichlet( v(1,1,1,1) )

# create a GTR rate matrix
Q := fnGTR(er,pi))");
	help_strings[string("fnGTR")][string("name")] = string(R"(fnGTR)");
	help_references[string("fnGTR")].push_back(RbHelpReference(R"(Tavare, S (1986). "Some Probabilistic and Statistical Problems in the Analysis of DNA Sequences".  Lectures on Mathematics in the Life Sciences. 17:57-86)",R"()",R"(http://www.damtp.cam.ac.uk/user/st321/CV_&_Publications_files/STpapers-pdf/T86.pdf)"));
	help_arrays[string("fnGTR")][string("see_also")].push_back(string(R"(fnJC)"));
	help_arrays[string("fnGTR")][string("see_also")].push_back(string(R"(fnK80)"));
	help_arrays[string("fnGTR")][string("see_also")].push_back(string(R"(fnK81)"));
	help_arrays[string("fnGTR")][string("see_also")].push_back(string(R"(fnT92)"));
	help_arrays[string("fnGTR")][string("see_also")].push_back(string(R"(fnHKY)"));
	help_arrays[string("fnGTR")][string("see_also")].push_back(string(R"(fnTrN)"));
	help_strings[string("fnGTR")][string("title")] = string(R"(The General Time-Reversible rate matrix)");
	help_arrays[string("fnGammaASRV")][string("authors")].push_back(string(R"(Benjamin Redelings)"));
	help_strings[string("fnGammaASRV")][string("description")] = string(R"(Add Gamma-distributed across-site rate variation (ASRV) to a site model.)");
	help_strings[string("fnGammaASRV")][string("details")] = string(R"(Each site evolves according to the specified site model, but at an unknown rate
that is Gamma distributed. If the site model parameter is a mixture model with
m components, this function will return a mixture with m*n components.

The continuous Gamma distribution is approximated with a mixture distribution
over n discrete rates, each with probability 1/n.  The Gamma distribution is
constrained to have a mean of 1, so as not to change the  branch lengths.
It therefore has only a single parameter alpha -- the shape parameter.
        - As alpha approaches infinity, all rates across sites become equal (rate variation goes to 0).
        - If alpha = 1, then the rate is exponentially distributed.  Rate variation is substantial.
        - As alpha approaches zero, many sites have rate 0, and many sites have a high rate.

RateMatrix and RateGenerator site model parameters will automatically be converted to a
SiteMixtureModel with a single component.)");
	help_strings[string("fnGammaASRV")][string("example")] = string(R"(# fnGammaASRV( ) constructs a mixture model that represents both the underlying
#   rate matrix and Gamma-distributed rate variation.
for (i in 1:10) { taxa[i] = taxon("T"+i) }
psi ~ dnBDP(lambda=1, rootAge=1, taxa=taxa)
alpha ~ dnExp(1/10)
er ~ dnDirichlet( [1,1,1,1,1,1] )
pi ~ dnDirichlet( [1,1,1,1] )
M := fnGammaASRV( fnGTR(er, pi), alpha, 4)
seq ~ dnPhyloCTMC(psi, M, type="DNA",nSites=10)

# As an alternative approach, models can be built up iteratively using pipes.
M := fnGTR(er,pi) |> fnGammaASRV(alpha, 4)

M := fnGTR(er,pi) |> fnGammaASRV(alpha, 4) |> fnInvASRV(p_inv)  # This has 5 (4+1) components - faster.
M := fnGTR(er,pi) |> fnInvASRV(p_inv) |> fnGammaASRV(alpha, 4)  # This has 8 (2*4) components - slower.

# The site model parameter can be a mixture model
weights ~ dnDirichlet([1,1])
pi1 ~ dnDirichlet( [1,1,1,1,1,1 ] )
pi2 ~ dnDirichlet( [1,1,1,1,1,1 ] )
M := fnMixtureASRV([fnGTR(er,pi1),fnGTR(er,pi2)],weights) |> fnGammaASRV(alpha) |> fnInvASRV(p_inv))");
	help_strings[string("fnGammaASRV")][string("name")] = string(R"(fnGammaASRV)");
	help_references[string("fnGammaASRV")].push_back(RbHelpReference(R"(Yang, Z. (1994) Maximum likelihood phylogenetic estimation from DNA sequences with variable rates over sites: approximate methods)",R"(https://doi.org/10.1007/BF00160154 )",R"()"));
	help_arrays[string("fnGammaASRV")][string("see_also")].push_back(string(R"(fnUnitMixture)"));
	help_arrays[string("fnGammaASRV")][string("see_also")].push_back(string(R"(fnInvASRV)"));
	help_arrays[string("fnGammaASRV")][string("see_also")].push_back(string(R"(fnScale)"));
	help_arrays[string("fnGammaASRV")][string("see_also")].push_back(string(R"(fnMixtureASRV)"));
	help_arrays[string("fnGammaASRV")][string("see_also")].push_back(string(R"(fnDiscretizeGamma)"));
	help_strings[string("fnGammaASRV")][string("title")] = string(R"(fnGammaASRV)");
	help_strings[string("fnGeographicalDistance")][string("name")] = string(R"(fnGeographicalDistance)");
	help_strings[string("fnHKY")][string("description")] = string(R"(DNA evolution model proposed in Hasegawa, Kishino, and Yano (1985).)");
	help_strings[string("fnHKY")][string("details")] = string(R"(In this model, nucleotides have different stationary frequencies, and transition and transversion rates are allowed to be different. Its first parameter, kappa, codes for the ratio between the rate of transitions and transversions. Its second parameter, baseFrequencies, codes for the frequencies of each nucleotide.

The HKY rate matrix elements will be of the form:
    Q[i, j] = c * kappa * baseFrequencies[j], if i<->j is a transition 
            = c * baseFrequencies[j], if i<->j is a transversion

where c is a constant needed to normalize the average rate to 1.)");
	help_strings[string("fnHKY")][string("example")] = string(R"(# the ratio between rates of transitions and transversions
kappa ~ dnLognormal(0,1)

# the base frequencies
pi ~ dnDirichlet( v(1,1,1,1) )

# create an HKY rate matrix
Q := fnHKY(kappa,pi))");
	help_strings[string("fnHKY")][string("name")] = string(R"(fnHKY)");
	help_references[string("fnHKY")].push_back(RbHelpReference(R"(Hasegawa, M. et al. (1985). "Dating of the human-ape splitting by a molecular clock of mitochondrial DNA". Journal of molecular evolution. 22(2):160-174.)",R"(https://doi.org/10.1007/BF02101694)",R"(https://link.springer.com/article/10.1007%2FBF02101694 )"));
	help_arrays[string("fnHKY")][string("see_also")].push_back(string(R"(fnJC)"));
	help_arrays[string("fnHKY")][string("see_also")].push_back(string(R"(fnK80)"));
	help_arrays[string("fnHKY")][string("see_also")].push_back(string(R"(fnK81)"));
	help_arrays[string("fnHKY")][string("see_also")].push_back(string(R"(fnF81)"));
	help_arrays[string("fnHKY")][string("see_also")].push_back(string(R"(fnT92)"));
	help_arrays[string("fnHKY")][string("see_also")].push_back(string(R"(fnTrN)"));
	help_arrays[string("fnHKY")][string("see_also")].push_back(string(R"(fnGTR)"));
	help_strings[string("fnHKY")][string("title")] = string(R"(The Hasegawa-Kishino-Yano (1985) nucleotide rate matrix)");
	help_strings[string("fnHiddenStateRateMatrix")][string("name")] = string(R"(fnHiddenStateRateMatrix)");
	help_strings[string("fnHostSwitchRateModifier")][string("name")] = string(R"(fnHostSwitchRateModifier)");
	help_strings[string("fnInfiniteSites")][string("name")] = string(R"(fnInfiniteSites)");
	help_arrays[string("fnInvASRV")][string("authors")].push_back(string(R"(Benjamin Redelings)"));
	help_strings[string("fnInvASRV")][string("description")] = string(R"(Add an invariable-sites component to a site model.)");
	help_strings[string("fnInvASRV")][string("details")] = string(R"(This model specifies that some fraction pInv of sites are invariable.
If the site model parameter is a mixture model with m components, this function will return a model with
m+1 components.)");
	help_strings[string("fnInvASRV")][string("example")] = string(R"(# fnInvASRV( ) creates a mixture model by adding invariant sites to an underlying site model.
for (i in 1:10) { taxa[i] = taxon("T"+i) }
psi ~ dnBDP(lambda=1, rootAge=1, taxa=taxa)
p_inv ~ dnUniform(0,1)
M := fnInvASRV( fnJC(4), p_inv)
seq ~ dnPhyloCTMC(psi, M, type="DNA", nSites=10)

# As an alternative approach, models can be built up iteratively using pipes.
M := fnJC(4) |> fnInv(p_inv)

M := fnJC(4) |> fnGammaASRV(alpha, 4) |> fnInvASRV(p_inv)  # This has 5 (4+1) components - faster.
M := fnJC(4) |> fnInvASRV(p_inv) |> fnGammaASRV(alpha, 4)  # This has 8 (4*2) components - slower.

# Not recommended -- illustration only.  3 components.
M := fnJC(4) |> fnInv(p1) |> fnInv(p2) # Fraction of invariable sites is p2 + (1-p2)*p1)");
	help_strings[string("fnInvASRV")][string("name")] = string(R"(fnInvASRV)");
	help_arrays[string("fnInvASRV")][string("see_also")].push_back(string(R"(fnUnitMixture)"));
	help_arrays[string("fnInvASRV")][string("see_also")].push_back(string(R"(fnGammaASRV)"));
	help_arrays[string("fnInvASRV")][string("see_also")].push_back(string(R"(fnMixtureASRV)"));
	help_arrays[string("fnInvASRV")][string("see_also")].push_back(string(R"(fnScale)"));
	help_strings[string("fnInvASRV")][string("title")] = string(R"(fnInvASRV)");
	help_strings[string("fnJC")][string("name")] = string(R"(fnJC)");
	help_strings[string("fnJones")][string("name")] = string(R"(fnJones)");
	help_strings[string("fnK80")][string("description")] = string(R"(DNA evolution model proposed in Kimura (1980).)");
	help_strings[string("fnK80")][string("details")] = string(R"(In this model, all nucleotides have an equal stationary frequency, and transition and transversion rates are allowed to be different. Its only parameter, kappa, codes for the ratio between the rate of transitions and transversions.

The K80 rate matrix elements will be of the form:
    Q[i, j] = c * kappa, if i<->j is a transition
            = c, if i<->j is a transversion

where c is a constant needed to normalize the average rate to 1.)");
	help_strings[string("fnK80")][string("example")] = string(R"(# the ratio between rates of transitions and transversions
kappa ~ dnExp(0.5)

# create a K80 rate matrix
Q := fnK80(kappa))");
	help_strings[string("fnK80")][string("name")] = string(R"(fnK80)");
	help_references[string("fnK80")].push_back(RbHelpReference(R"(Kimura M (1980). "A simple method for estimating evolutionary rates of base substitutions through comparative studies of nucleotide sequences". Journal of Molecular Evolution. 16:111–20.)",R"(https://doi.org/10.1007/BF01731581)",R"(https://link.springer.com/article/10.1007/BF01731581 )"));
	help_arrays[string("fnK80")][string("see_also")].push_back(string(R"(fnJC)"));
	help_arrays[string("fnK80")][string("see_also")].push_back(string(R"(fnF81)"));
	help_arrays[string("fnK80")][string("see_also")].push_back(string(R"(fnK81)"));
	help_arrays[string("fnK80")][string("see_also")].push_back(string(R"(fnT92)"));
	help_arrays[string("fnK80")][string("see_also")].push_back(string(R"(fnHKY)"));
	help_arrays[string("fnK80")][string("see_also")].push_back(string(R"(fnGTR)"));
	help_strings[string("fnK80")][string("title")] = string(R"(The Kimura (1980) nucleotide rate matrix)");
	help_strings[string("fnK81")][string("description")] = string(R"(DNA evolution model proposed in Kimura (1981).)");
	help_strings[string("fnK81")][string("details")] = string(R"(In this model, transition and transversion rates are allowed to be different, and transversion rates for A <-> C, G <-> T and A <-> T, C <-> G transversions are different as well. The first argument, kappa1, defines the ratio between the rate of transitions and the rate of A <-> C, G <-> T transversions. The second argument, kappa2, defines the ratio between the rate of A <-> T, C <-> G transversions and the rate of A <-> C, G <-> T transversions. The third argument, baseFrequencies, defines the stationary frequencies of nucleotide bases. Note that the original Kimura (1981) model assumed equal base frequencies, so this function is more general (if ran without a baseFrequencies argument, however, this is equivalent to K81, since the default is all frequencies equal). 

The K81 rate matrix elements will be of the form:
    Q[i, j] = c, if i<->j is an A<->C/G<->T transversion
            = c * kappa1, if i<->j is a transition
            = c * kappa2, if i<->j is an A<->T/C<->G transversion

where c is a constant needed to normalize the average rate to 1. If using the baseFrequencies parameter, those elements are multiplied by baseFrequencies[j].)");
	help_strings[string("fnK81")][string("example")] = string(R"(# the ratio between rates of transitions and A<->C/G<->T transversions
kappa1 ~ dnExp(0.5)

# the ratio between rates of A<->T/C<->G and A<->C/G<->T transversions
kappa2 ~ dnExp(0.5)

# create a K81 rate matrix
Q := fnK81(kappa1, kappa2)

# base frequencies
baseFrequencies ~ dnDirichlet(v(1,1,1,1))

# K81 rate matrix with non-equal base frequencies
Q := fnK81(kappa1, kappa2, baseFrequencies))");
	help_strings[string("fnK81")][string("name")] = string(R"(fnK81

# title
The Kimura (1981) nucleotide rate matrix)");
	help_references[string("fnK81")].push_back(RbHelpReference(R"(Kimura M (1981). "Estimation of evolutionary distances between homologous nucleotide sequences". Proceedings of the National Academy of Sciences of the United States of America. 78:454–8.)",R"(https://doi.org/10.1073/pnas.78.1.454)",R"(https://www.pnas.org/doi/abs/10.1073/pnas.78.1.454 )"));
	help_arrays[string("fnK81")][string("see_also")].push_back(string(R"(fnJC)"));
	help_arrays[string("fnK81")][string("see_also")].push_back(string(R"(fnK80)"));
	help_arrays[string("fnK81")][string("see_also")].push_back(string(R"(fnF80)"));
	help_arrays[string("fnK81")][string("see_also")].push_back(string(R"(fnT92)"));
	help_arrays[string("fnK81")][string("see_also")].push_back(string(R"(fnHKY)"));
	help_arrays[string("fnK81")][string("see_also")].push_back(string(R"(fnTrN)"));
	help_arrays[string("fnK81")][string("see_also")].push_back(string(R"(fnGTR)"));
	help_strings[string("fnLG")][string("name")] = string(R"(fnLG)");
	help_strings[string("fnLnProbability")][string("name")] = string(R"(fnLnProbability)");
	help_arrays[string("fnMinBLTimeScaling")][string("authors")].push_back(string(R"(David Černý)"));
	help_arrays[string("fnMinBLTimeScaling")][string("authors")].push_back(string(R"(Laura Mulvey)"));
	help_strings[string("fnMinBLTimeScaling")][string("description")] = string(R"(Time-scales an undated tree based on a vector of tip ages using the minimum
branch length ("MBL") approach (Laurin 2004; Bapst 2014).)");
	help_strings[string("fnMinBLTimeScaling")][string("details")] = string(R"(The age of each internal node is based on the age of the oldest tip descended
from it. However, if t0 is the oldest tip descended from a given node (denoted
x) and also the oldest tip descended from that node's parent (denoted y), then
setting t(x) = t(y) = t(t0) would produce zero-length branches y->x and x->t0.
We avoid this by requiring every branch to be no shorter than some user-supplied
constant. This has the effect of shifting node ages deeper into the past.

Conceptually, the undated tree would usually correspond either to a bare
topology (a tree without branch lengths) or a tree with branch lengths in units
of expected change; in practice, both `BranchLengthTree` and `TimeTree` arguments
are accepted. In this implementation of the MBL approach, both terminal and
internal branches are required to be greater than or equal to the specified
minimum. If there is uncertainty associated with the age of a given tip,
the midpoint of the uncertainty range is used for time-scaling.

The algorithm is not stochastic (i.e., it always returns the same time-scaled
tree for a given input), and is primarily intended to generate a plausible
starting tree for MCMC analyses.)");
	help_strings[string("fnMinBLTimeScaling")][string("example")] = string(R"(# Read in an undated tree
undated_tree <- readTrees("undated.nex")[1]

# Read tip age data from a file
taxa <- readTaxonData("tipages.tsv")

# Time-scale using a minimum branch length of 3 Myr
dated_tree <- fnMinBLTimeScaling(undated_tree, taxa, 3.0)

print(undated_tree) # The original tree remains unchanged
print(dated_tree)   # A new, dated tree has been returned)");
	help_strings[string("fnMinBLTimeScaling")][string("name")] = string(R"(fnMinBLTimeScaling)");
	help_references[string("fnMinBLTimeScaling")].push_back(RbHelpReference(R"(Bapst DW (2014). Assessing the effect of time-scaling methods on phylogeny-based analyses in the fossil record. Paleobiology, 40(3):331-351.)",R"(10.1666/13033)",R"()"));
	help_references[string("fnMinBLTimeScaling")].push_back(RbHelpReference(R"(Laurin M (2004). The evolution of body size, Cope's rule and the origin of amniotes. Systematic Biology, 53(4):594-622.)",R"(10.1080/10635150490445706)",R"()"));
	help_arrays[string("fnMinBLTimeScaling")][string("see_also")].push_back(string(R"(simStartingTree)"));
	help_arrays[string("fnMixtureASRV")][string("authors")].push_back(string(R"(Benjamin Redelings)"));
	help_strings[string("fnMixtureASRV")][string("description")] = string(R"(Constructs a mixture model from a collection of site models.)");
	help_strings[string("fnMixtureASRV")][string("details")] = string(R"(Each site will evolve according to one of the input site models, which may also
be mixture models.  The probability that each site follows a particular site model
is specified by the fractions parameter.

The number of components in the resulting mixture model is the sum of the number
of components of the input mixture models.

If the fractions parameter is missing, then each of the given models is given equal
weight.)");
	help_strings[string("fnMixtureASRV")][string("example")] = string(R"(# Two components with different frequencies
for (i in 1:10) { taxa[i] = taxon("T"+i) }
psi ~ dnBDP(lambda=1, rootAge=1, taxa=taxa)
pi1 ~ dnDirichlet([1,1,1,1])
pi2 ~ dnDirichlet([1,1,1,1])
weights ~ dnDirichlet([1,1])
M := fnMixtureASRV([fnF81(pi1),fnF81(pi2)],weights)
seq ~ dnPhyloCTMC(psi, M, type="DNA", nSites=10)

# A weight of 1/2 on each model because the weights are missing.
M := fnMixtureASRV([fnF81(pi1),fnF81(pi2)])

# Adding rate variation to the frequency-variation model.
M := fnMixtureASRV([fnF81(pi1),fnF81(pi2)],weights) |> fnGammaASRV(alpha) |> fnInvASRV(p_inv))");
	help_strings[string("fnMixtureASRV")][string("name")] = string(R"(fnMixtureASRV)");
	help_arrays[string("fnMixtureASRV")][string("see_also")].push_back(string(R"(fnUnitMixture)"));
	help_arrays[string("fnMixtureASRV")][string("see_also")].push_back(string(R"(fnGammaASRV)"));
	help_arrays[string("fnMixtureASRV")][string("see_also")].push_back(string(R"(fnInvASRV)"));
	help_arrays[string("fnMixtureASRV")][string("see_also")].push_back(string(R"(fnScale)"));
	help_strings[string("fnMixtureASRV")][string("title")] = string(R"(fnMixtureASRV)");
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
	help_arrays[string("fnScale")][string("authors")].push_back(string(R"(Benjamin Redelings)"));
	help_strings[string("fnScale")][string("description")] = string(R"(Scale a vector of SiteMixtureModels)");
	help_strings[string("fnScale")][string("details")] = string(R"(This function has two forms.  The first form takes a SiteMixtureModel `model` and scales it by
a rate `rate`.  This form returns SiteMixtureModel.

The second form takes SiteMixtureModel[] `models` and RealPos[] `rates`, and scales `models[i]`
by `rates[i]`.  This form returns SiteMixtureModel[].

As a shortcut, if the second argument `rates` is a vector but the first element `model` is not,
then the first argument will be automatically replaced with a vector of SiteMixtureModels of the
same length as `rates`, where each element is identical to `model`.)");
	help_strings[string("fnScale")][string("example")] = string(R"(Q = fnJC(4)                    # The rate of Q is 1

# Operating on SiteMixtureModel
Q2 = fnScale(Q,2)              # The rate of Q2 is 2

# Operating on SiteMixtureModel[]
Qs = fnScale([Q,Q],[1,2])      # Qs[1] and Qs[2] have rates 1 and 2
Qs = fnScale(Q,    [1,2])      # An abbreviation for the above.

# We can build up models iteratively using pipes
Qs = Q |> fnScale([1,2])       # A shorter abbreviation.

# A JC+LogNormal[4] ASRV model
site_rates := dnLognormal(0,lsigma) |> fnDiscretizeDistribution(4)
MM := fnJC(4) |> fnScale(site_rates) |> fnMixtureASRV()
M := fnScale(MM, 1/MM.rate())

# A FreeRates[5] ASRV model
rates ~ dnDirichlet( [1,1,1,1,1] )
weights ~ dnDirichlet( [2,2,2,2,2] )
MM := fnJC(4) |> fnScale(rates) |> fnMixtureASRV(weights)
M := fnScale(MM, 1/MM.rate()))");
	help_strings[string("fnScale")][string("name")] = string(R"(fnScale)");
	help_arrays[string("fnScale")][string("see_also")].push_back(string(R"(fnUnitMixture)"));
	help_arrays[string("fnScale")][string("see_also")].push_back(string(R"(fnInvASRV)"));
	help_arrays[string("fnScale")][string("see_also")].push_back(string(R"(fnMixtureASRV)"));
	help_strings[string("fnScale")][string("title")] = string(R"(fnScale)");
	help_strings[string("fnSegregatingSites")][string("name")] = string(R"(fnSegregatingSites)");
	help_strings[string("fnShortestDistance")][string("name")] = string(R"(fnShortestDistance)");
	help_strings[string("fnSiteRateModifier")][string("name")] = string(R"(fnSiteRateModifier)");
	help_arrays[string("fnSmoothTimeLine")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("fnSmoothTimeLine")][string("description")] = string(R"(Function to create a smooth timeline where all values after a maximum time are constant, i.e., equal to the previous interval, to avoid crazy looking plots from the prior.)");
	help_strings[string("fnSmoothTimeLine")][string("details")] = string(R"(Thus function takes a vector of values and a matching vector of times and a maximum time. Then, it constructs a smooth timeline by using all values before the maximum, and replacing all values after the maximum with the last value before the maximum. Thus, the timeline is smooth after the maximum.)");
	help_strings[string("fnSmoothTimeLine")][string("name")] = string(R"(fnSmoothTimeLine)");
	help_strings[string("fnSmoothTimeLine")][string("title")] = string(R"(Create a smooth timeline)");
	help_strings[string("fnStateCountRateModifier")][string("name")] = string(R"(fnStateCountRateModifier)");
	help_strings[string("fnStirling")][string("name")] = string(R"(fnStirling)");
	help_strings[string("fnStitchTree")][string("name")] = string(R"(fnStitchTree)");
	help_strings[string("fnT92")][string("description")] = string(R"(DNA evolution model proposed in Tamura (1992).)");
	help_strings[string("fnT92")][string("details")] = string(R"(In this model, A and T have an equal stationary frequency, with G and C frequencies distinct, and transition and transversion rates are allowed to be different. Its first parameter, kappa, codes for the ratio between the rate of transitions and transversions. Its second parameter, gc, codes for the compound frequency of G and C nucleotides.

The T92 rate matrix elements will be of the form:
    Q[i, j] = c * kappa * gc / 2, if i<->j is a transition and j is C or G
            = c * gc / 2, if i<->j is a transversion and j is C or G
            = c * kappa * (1 - gc) / 2, if i<->j is a transition and j is A or T
            = c * (1 - gc) / 2, if i<->j is a transversion and j is A or T

where c is a constant needed to normalize the average rate to 1.)");
	help_strings[string("fnT92")][string("example")] = string(R"(# the ratio between rates of transitions and transversions
kappa ~ dnExp(0.5)

# the frequency of G and C nucleotides
gc ~ dnUnif(0, 1)

# create a T92 rate matrix
Q := fnT92(kappa, gc))");
	help_strings[string("fnT92")][string("name")] = string(R"(fnT92)");
	help_references[string("fnT92")].push_back(RbHelpReference(R"(Tamura K (1992). "Estimation of the number of nucleotide substitutions when there are strong transition-transversion and G+C-content biases". Molecular Biology and Evolution. 9:678–87.)",R"(https://doi.org/10.1093/oxfordjournals.molbev.a040752)",R"(https://academic.oup.com/mbe/article/9/4/678/1254082 )"));
	help_arrays[string("fnT92")][string("see_also")].push_back(string(R"(fnJC)"));
	help_arrays[string("fnT92")][string("see_also")].push_back(string(R"(fnF81)"));
	help_arrays[string("fnT92")][string("see_also")].push_back(string(R"(fnK80)"));
	help_arrays[string("fnT92")][string("see_also")].push_back(string(R"(fnK81)"));
	help_arrays[string("fnT92")][string("see_also")].push_back(string(R"(fnHKY)"));
	help_arrays[string("fnT92")][string("see_also")].push_back(string(R"(fnTrN)"));
	help_arrays[string("fnT92")][string("see_also")].push_back(string(R"(fnGTR)"));
	help_strings[string("fnT92")][string("title")] = string(R"(The Tamura (1992) nucleotide rate matrix)");
	help_strings[string("fnTIM")][string("name")] = string(R"(fnTIM)");
	help_strings[string("fnTVM")][string("name")] = string(R"(fnTVM)");
	help_strings[string("fnTajimasD")][string("name")] = string(R"(fnTajimasD)");
	help_strings[string("fnTajimasPi")][string("name")] = string(R"(fnTajimasPi)");
	help_strings[string("fnTrN")][string("description")] = string(R"(DNA evolution model proposed in Tamura & Nei (1993).)");
	help_strings[string("fnTrN")][string("details")] = string(R"(In this model, nucleotide base frequencies are different, and the two transition rates (A <-> G and C<->T) can be different to each other, and to the transversion rate. The first argument, kappa1, defines the ratio between the rate of A <-> G (i.e. purine) transitions to transversions. The second argument, kappa2, defines the ratio between the rate of C <-> T (i.e. pyrimidine) transitions to transversions. The third argument, baseFrequencies, defines the stationary frequencies of nucleotide bases. 

The TrN rate matrix elements are of the form:
    Q[i, j] = c * kappa1 * baseFrequencies[j], if i<->j is A<->G
            = c * kappa2 * baseFrequencies[j], if i<->j is C<->T
            = c * baseFrequencies[j], otherwise

where c is a constant needed to normalize the average rate to 1)");
	help_strings[string("fnTrN")][string("example")] = string(R"(# A <-> G transition rate
kappaAG ~ dnLognormal(0,1)

# C <-> T transition rate
kappaCT ~ dnLognormal(0,1)

# nucleotide base frequencies
pi ~ dnDirichlet( v(1,1,1,1) )

# create a TrN rate matrix
Q := fnTrN(kappaAG, kappaCT, ,pi))");
	help_strings[string("fnTrN")][string("name")] = string(R"(fnTrN)");
	help_references[string("fnTrN")].push_back(RbHelpReference(R"(Tamura, K. and M. Nei (1993). "Estimation of the number of nucleotide substitutions in the control region of mitochondrial DNA in humans and chimpanzees". Molecular biology and evolution. 10(3):512-526.)",R"(https://doi.org/10.1093/oxfordjournals.molbev.a040023)",R"(https://academic.oup.com/mbe/article/10/3/512/1016366 )"));
	help_arrays[string("fnTrN")][string("see_also")].push_back(string(R"(fnJC)"));
	help_arrays[string("fnTrN")][string("see_also")].push_back(string(R"(fnK80)"));
	help_arrays[string("fnTrN")][string("see_also")].push_back(string(R"(fnK81)"));
	help_arrays[string("fnTrN")][string("see_also")].push_back(string(R"(fnT92)"));
	help_arrays[string("fnTrN")][string("see_also")].push_back(string(R"(fnHKY)"));
	help_arrays[string("fnTrN")][string("see_also")].push_back(string(R"(fnGTR)"));
	help_strings[string("fnTrN")][string("title")] = string(R"(The Tamura-Nei (1993) nucleotide rate matrix)");
	help_strings[string("fnTreeAssembly")][string("name")] = string(R"(fnTreeAssembly)");
	help_strings[string("fnTreePairwiseDistances")][string("name")] = string(R"(fnTreePairwiseDistances)");
	help_strings[string("fnTreePairwiseNodalDistances")][string("name")] = string(R"(fnTreePairwiseNodalDistances)");
	help_strings[string("fnTreeScale")][string("name")] = string(R"(fnTreeScale)");
	help_strings[string("fnUnitMixture")][string("description")] = string(R"(Create a SiteMixtureModel from a RateMatrix or RateGenerator)");
	help_strings[string("fnUnitMixture")][string("details")] = string(R"(This function creates a SiteMixtureModel with one component by specifying the
rate and root frequencies for a RateGenerator.  The rate defaults to 1, leaving
the underlying model unchanged.

If the site model parameter is a RateMatrix, the root frequencies default to the
equilibrium frequencies of the RateMatrix.  However, a RateGenerator might not have
equilibrium frequencies, in which case the root frequencies must be specified explicitly.

In many cases it is not necessary to explicitly call fnUnitMixture(), RevBayes can
automatically convert a RateMatrix to a SiteMixtureModel.)");
	help_strings[string("fnUnitMixture")][string("example")] = string(R"(M := fnUnitMixture( fnJC(4) )
M := fnJC(4) |> fnUnitMixture()  # nested functions can be expressed using pipes.

# Explicit conversion to SiteMixtureModel
M := fnGTR(er,pi) |> fnUnitMixture() |> fnGammaASRV(alpha) |> fnInvASRV(p_inv)
# Implicit conversion to SiteMixtureModel
M := fnGTR(er,pi) |> fnGammaASRV(alpha) |> fnInvASRV(p_inv)

# Specifying the root frequencies
M := fnDECRateMatrix(dr,er,"Include") |> fnUnitMixture(rootFrequencies=simplex(rep(1,n_states))))");
	help_strings[string("fnUnitMixture")][string("name")] = string(R"(fnUnitMixture)");
	help_strings[string("fnUnitMixture")][string("title")] = string(R"(fnUnitMixture)");
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
	help_strings[string("getOption")][string("details")] = string(R"(Runtime options are used to personalize RevBayes and are stored on the local machine. See `setOption` for the list of available keys and their associated values.)");
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
	help_strings[string("logistic")][string("description")] = string(R"(Compute the logistic function)");
	help_strings[string("logistic")][string("details")] = string(R"(The function is defined as

        logistic(x) = 1/(1 + exp(-x))

                    = exp(x)/(1 + exp(x))

This function takes a real number to a probability.
It is the inverse of the logit function.)");
	help_strings[string("logistic")][string("example")] = string(R"(x ~ dnNormal(0,1)
p := logistic(x))");
	help_strings[string("logistic")][string("name")] = string(R"(logistic)");
	help_arrays[string("logistic")][string("see_also")].push_back(string(R"(logit)"));
	help_strings[string("logistic")][string("title")] = string(R"(The logistic function)");
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
	help_strings[string("mcmc")][string("details")] = string(R"(The MCMC analysis object produced by a call to this function keeps copies of the model and the associated moves and monitors. The MCMC analysis object is used to run Markov chain Monte Carlo (MCMC) simulation on the model, using the provided moves, to obtain a sample of the posterior probability distribution. During the analysis, the monitors are responsible for sampling model parameters of interest.

The `mcmc.run()` method begins or continues an MCMC analysis. 

During each iteration of an analysis, moves are selected from those listed in the `moves` parameter.  With the default `moveschedule = "random"`, or `moveschedule = "sequential"`, moves will be attempted, on average, `weight` times per iteration.  If `moveschedule = "single"`, RevBayes will attempt exactly one move per iteration, corresponding to the behavior of software like BEAST or MrBayes. See Höhna et al. (2017) for details.

The run will continue for `generations` iterations, or until a stopping rule is triggered: perhaps once the run has attained convergence, or after a certain amount of time has passed.  The run will be terminated once *all* convergence rules ([`srGelmanRubin()`], [`srGeweke()`], [`srMinESS()`], [`srStationarity()`]) in its `StoppingRule[]` argument have been fulfilled; or once *any* threshold rules ([`srMaxTime()`], [`srMaxIteration()`]) are met.

The parameters `checkpointFile` and `checkpointInterval` generate snapshots of the current state of the MCMC run from which the run can be continued if interrupted using the `mcmc.initializeFromCheckpoint()` method.

The `mcmc.initializeFromCheckpoint()` method allows an analysis to be continued from a checkpoint file. New generations will be appended to existing monitor files.)");
	help_strings[string("mcmc")][string("example")] = string(R"(# Create a simple model (unclamped)
a ~ dnExponential(1)
mymodel = model(a)

# Create a move vector and a monitor vector
moves[1] = mvScale( a, lambda = 1.0, weight = 1.0 )
monitors[1] = mnFile( a, filename = "output/out.log" )

# Create an mcmc object
mymcmcObject = mcmc( mymodel, monitors, moves )

# Run a short analysis
mymcmcObject.burnin( generations = 400, tuningInterval = 100 )
mymcmcObject.run( generations = 400, checkpointFile = "output/out.ckp", checkpointInterval = 100 )

# print the summary of the operators (now tuned)
mymcmcObject.operatorSummary()

# Resume analysis from the checkpoint file
mymcmcObject.initializeFromCheckpoint( "output/out.ckp" )

# Conduct an additional 400 generations
mymcmcObject.run( generations = 400 )

# Stopping rules are defined on the total number of generations
# This command will have no effect, as 400 generations have already been performed.
mymcmcObject.run( rules = [ srMaxIteration(400) ] ))");
	help_strings[string("mcmc")][string("name")] = string(R"(mcmc)");
	help_references[string("mcmc")].push_back(RbHelpReference(R"(Metropolis N, AW Rosenbluth, MN Rosenbluth, AH Teller, E Teller (1953). Equation of state calculations by fast computing machines. Journal of Chemical Physics, 21:1087-1092.)",R"(10.1063/1.1699114)",R"()"));
	help_references[string("mcmc")].push_back(RbHelpReference(R"(Hastings WK (1970) Monte Carlo sampling methods using Markov chains and their applications. Biometrika, 57:97-109.)",R"(10.2307/2334940)",R"()"));
	help_references[string("mcmc")].push_back(RbHelpReference(R"(Höhna S, Landis MJ, Heath TA (2017).  Phylogenetic inference using `RevBayes`. Current Protocols in Bioinformatics.)",R"(10.1002/cpbi.22)",R"()"));
	help_arrays[string("mcmc")][string("see_also")].push_back(string(R"(mcmcmc)"));
	help_strings[string("mcmc")][string("title")] = string(R"(MCMC analysis object)");
	help_arrays[string("mcmcmc")][string("authors")].push_back(string(R"(Michael Landis)"));
	help_arrays[string("mcmcmc")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("mcmcmc")][string("description")] = string(R"(The Mcmcmc analysis object keeps a model and the associated moves and monitors. The object is used to run Metropolis Couped Markov chain Monte Carlo (Mcmcmc) simulation on the model, using the provided moves, to obtain a sample of the posterior probability distribution. During the analysis, the monitors are responsible for sampling model parameters of interest.)");
	help_strings[string("mcmcmc")][string("details")] = string(R"(The Mcmcmc analysis object produced by a call to this function keeps copies of the model and the associated moves and monitors. The Mcmcmc analysis object is used to run Markov chain Monte Carlo (Mcmcmc) simulation on the model, using the provided moves, to obtain a sample of the posterior probability distribution. During the analysis, the monitors are responsible for sampling model parameters of interest.

An MCMCMC analysis is initiated using the `mcmcmc.run()` method.  
The `StoppingRule[]` argument provides a mechanism to automatically terminate a run once a set of rules are met: perhaps once the run has attained convergence, or after a certain amount of time has passed.  The run will be terminated once *all* convergence rules ([`srGelmanRubin()`], [`srGeweke()`], [`srMinESS()`], [`srStationarity()`]) have been fulfilled; or once *any* threshold rules ([`srMaxTime()`], [`srMaxIteration()`]) are met.
The parameters `checkpointFile` and `checkpointInterval` generate snapshots of the current state of the MCMCMC run from which the run can be continued if interrupted using the `mcmc.initializeFromCheckpoint()` method. An example is given on the documentation page for [`mcmc()`].)");
	help_strings[string("mcmcmc")][string("example")] = string(R"(# Create a simple model (unclamped)
a ~ dnExponential(1)
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
	help_strings[string("model")][string("description")] = string(R"(Creates a model object that can be graphed or subjected to Bayesian inference.)");
	help_strings[string("model")][string("details")] = string(R"(`model(x)` creates a model object by creating a copy of all elements and 
parameters that influence or are influenced by the likelihood of `x`.

Because `model` works with copies of objects, conducting an mcmc(mc) analysis
on a model object will not change the values of the objects in the RevBayes
workspace.

The model object can be modified to ignore specific data elements using the
method `ignoreData`.  Thus to run without the sequence data `phySeq` you
might specify:

   mymodel.ignoreData(phySeq)

Only clamped nodes can be ignored. To ignore all clamped nodes you can use
the method `ignoreAllData`:

   mymodel.ignoreAllData())");
	help_strings[string("model")][string("example")] = string(R"(# Create a simple model (unclamped)
a ~ dnExponential(1)
b ~ dnExponential(a)
mymodel = model(b) # equivalent to model(a) or model(a, b)

# Save a DOT visualization of the model to file
mymodel.graph("mymodel.dot")

# Create a move vector and a monitor vector
moves = [ mvScale( a, lambda = 1.0, weight = 1.0 ) ]
monitors = [ mnScreen(printgen = 10, a) ]

# Create an mcmc object
mymcmcObject = mcmc( mymodel, monitors, moves )

# Print value of a
print(a)

# Run a short analysis
mymcmcObject.run( generations = 100 )

print(a) # Value is unchanged in the workspace - only the copy is modified)");
	help_strings[string("model")][string("name")] = string(R"(model)");
	help_strings[string("model")][string("title")] = string(R"(Create a model object)");
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
	help_strings[string("mvBetaSimplex")][string("description")] = string(R"(The Beta Simplex move selects one element of the a vector and proposes a new value for it drawn from a Beta distribution.)");
	help_strings[string("mvBetaSimplex")][string("example")] = string(R"(A usage example can be found at https://revbayes.github.io/tutorials/chromo/#root)");
	help_strings[string("mvBetaSimplex")][string("name")] = string(R"(mvBetaSimplex)");
	help_arrays[string("mvBetaSimplex")][string("see_also")].push_back(string(R"(- mvDirichletSimplex)"));
	help_arrays[string("mvBetaSimplex")][string("see_also")].push_back(string(R"(- mvElementSwapSimplex)"));
	help_strings[string("mvBetaSimplex")][string("title")] = string(R"(Beta Simplex move)");
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
	help_arrays[string("mvDPPValueBetaSimplex")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("mvDPPValueBetaSimplex")][string("description")] = string(R"(Operates on draws from a Dirichlet process prior (DPP) on mixtures of [Simplex](https://revbayes.github.io/documentation/Simplex.html) distributions, i.e., distributions defined over vectors whose elements are positive and sum to 1.)");
	help_strings[string("mvDPPValueBetaSimplex")][string("details")] = string(R"(In Dirichlet process mixtures, the number of categories (= clusters) is not specified beforehand but inferred from the data, and can range anywhere from 1 to the total number of elements (= observations). The move takes the current number of categories and simultaneously updates the value of every category using the beta simplex move with a concentration parameter (alpha) of 10.)");
	help_strings[string("mvDPPValueBetaSimplex")][string("example")] = string(R"(# Here, we draw from a DP mixture for 3 elements, where every element
# is itself a 2-element simplex drawn from a flat Dirichlet distribution
x ~ dnDPP(1, dnDirichlet( [1, 1] ), 3)

# Next, we add the move. Note that without moves other than
# mvDPPValueBetaSimplex, only the values of the categories will be
# updated: the total number of categories and the assignment of elements
# to categories will be determined by the initial draw.
moves[1] = mvDPPValueBetaSimplex(x, weight=1)

monitors[1] = mnScreen(x, printgen=1)
mymodel = model(x)
mymcmc = mcmc(mymodel, monitors, moves)
mymcmc.run(generations=50))");
	help_strings[string("mvDPPValueBetaSimplex")][string("name")] = string(R"(mvDPPValueBetaSimplex)");
	help_arrays[string("mvDPPValueBetaSimplex")][string("see_also")].push_back(string(R"(- dnDPP)"));
	help_arrays[string("mvDPPValueBetaSimplex")][string("see_also")].push_back(string(R"(- mvBetaSimplex)"));
	help_arrays[string("mvDPPValueBetaSimplex")][string("see_also")].push_back(string(R"(- mvDPPValueScaling)"));
	help_arrays[string("mvDPPValueBetaSimplex")][string("see_also")].push_back(string(R"(- mvDPPValueSliding)"));
	help_strings[string("mvDPPValueBetaSimplex")][string("title")] = string(R"(Beta simplex move applied to individual categories of a Dirichlet process mixture)");
	help_arrays[string("mvDPPValueScaling")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("mvDPPValueScaling")][string("description")] = string(R"(Operates on draws from a Dirichlet process prior (DPP) on mixtures of [RealPos](https://revbayes.github.io/documentation/RealPos.html) distributions, i.e., distributions defined over non-negative real numbers.)");
	help_strings[string("mvDPPValueScaling")][string("details")] = string(R"(In Dirichlet process mixtures, the number of categories (= clusters) is not specified beforehand but inferred from the data, and can range anywhere from 1 to the total number of elements (= observations). The move takes the current number of categories and simultaneously updates the value of every category using the scaling move with a scaling factor (lambda) of 1.)");
	help_strings[string("mvDPPValueScaling")][string("example")] = string(R"(# Here, we draw from a DP mixture for 3 elements, where every element
# is a non-negative real number drawn from an exponential distribution
x ~ dnDPP(1, dnExp(1), 3)

# Next, we add the move. Note that without moves other than
# mvDPPValueScaling, only the values of the categories will be updated:
# the total number of categories and the assignment of elements to
# categories will be determined by the initial draw.
moves[1] = mvDPPValueScaling(x, weight=1)

monitors[1] = mnScreen(x, printgen=1)
mymodel = model(x)
mymcmc = mcmc(mymodel, monitors, moves)
mymcmc.run(generations=50))");
	help_strings[string("mvDPPValueScaling")][string("name")] = string(R"(mvDPPValueScaling)");
	help_arrays[string("mvDPPValueScaling")][string("see_also")].push_back(string(R"(- dnDPP)"));
	help_arrays[string("mvDPPValueScaling")][string("see_also")].push_back(string(R"(- mvScale)"));
	help_arrays[string("mvDPPValueScaling")][string("see_also")].push_back(string(R"(- mvDPPValueBetaSimplex)"));
	help_arrays[string("mvDPPValueScaling")][string("see_also")].push_back(string(R"(- mvDPPValueSliding)"));
	help_strings[string("mvDPPValueScaling")][string("title")] = string(R"(Scaling move applied to individual categories of a Dirichlet process mixture)");
	help_arrays[string("mvDPPValueSliding")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("mvDPPValueSliding")][string("description")] = string(R"(Operates on draws from a Dirichlet process prior (DPP) on mixtures of [Real](https://revbayes.github.io/documentation/Real.html) distributions, i.e., distributions defined over all real numbers.)");
	help_strings[string("mvDPPValueSliding")][string("details")] = string(R"(In Dirichlet process mixtures, the number of categories (= clusters) is not specified beforehand but inferred from the data, and can range anywhere from 1 to the total number of elements (= observations). The move takes the current number of categories and simultaneously updates the value of every category using the sliding move with a window size (delta) of 1.)");
	help_strings[string("mvDPPValueSliding")][string("example")] = string(R"(# Here, we draw from a DP mixture for 3 elements, where every element
# is a real number drawn from the standard normal distribution
x ~ dnDPP(1, dnNormal(0, 1), 3)

# Next, we add the move. Note that without moves other than
# mvDPPValueSliding, only the values of the categories will be updated:
# the total number of categories and the assignment of elements to
# categories will be determined by the initial draw.
moves[1] = mvDPPValueSliding(x, weight=1)

monitors[1] = mnScreen(x, printgen=1)
mymodel = model(x)
mymcmc = mcmc(mymodel, monitors, moves)
mymcmc.run(generations=50))");
	help_strings[string("mvDPPValueSliding")][string("name")] = string(R"(mvDPPValueSliding)");
	help_arrays[string("mvDPPValueSliding")][string("see_also")].push_back(string(R"(- dnDPP)"));
	help_arrays[string("mvDPPValueSliding")][string("see_also")].push_back(string(R"(- mvSlide)"));
	help_arrays[string("mvDPPValueSliding")][string("see_also")].push_back(string(R"(- mvDPPValueBetaSimplex)"));
	help_arrays[string("mvDPPValueSliding")][string("see_also")].push_back(string(R"(- mvDPPValueScaling)"));
	help_strings[string("mvDPPValueSliding")][string("title")] = string(R"(Sliding move applied to individual categories of a Dirichlet process mixture)");
	help_strings[string("mvDirichletSimplex")][string("description")] = string(R"(A Dirichlet-simplex proposal randomly changes some values of a [Simplex](https://revbayes.github.io/documentation/Simplex.html)
(a vector whose elements sum to 1). The other values change too because of renormalization.
 
First, some random indices are drawn.
Then, the proposal draws a new simplex `u ~ Dirichlet(val[index] * alpha)`, where alpha is the tuning parameter.
The new value is set to `u`.
The simplex is then renormalized.)");
	help_strings[string("mvDirichletSimplex")][string("example")] = string(R"(Usage examples can be found at https://revbayes.github.io/tutorials/morph_tree/V2.html and https://revbayes.github.io/tutorials/morph_ase/ase_free.html)");
	help_strings[string("mvDirichletSimplex")][string("name")] = string(R"(mvDirichletSimplex)");
	help_arrays[string("mvDirichletSimplex")][string("see_also")].push_back(string(R"(- mvBetaSimplex)"));
	help_arrays[string("mvDirichletSimplex")][string("see_also")].push_back(string(R"(- mvElementSwapSimplex)"));
	help_strings[string("mvDirichletSimplex")][string("title")] = string(R"(Dirichlet Simplex move)");
	help_strings[string("mvDiscreteEventCategoryRandomWalk")][string("name")] = string(R"(mvDiscreteEventCategoryRandomWalk)");
	help_strings[string("mvElementSwapSimplex")][string("description")] = string(R"(An Element Swap Simplex move selects two elements of a vector and exchanges their values.)");
	help_strings[string("mvElementSwapSimplex")][string("name")] = string(R"(mvElementSwapSimplex)");
	help_arrays[string("mvElementSwapSimplex")][string("see_also")].push_back(string(R"(mvBetaSimplex)"));
	help_arrays[string("mvElementSwapSimplex")][string("see_also")].push_back(string(R"(mvDirichletSimplex)"));
	help_strings[string("mvElementSwapSimplex")][string("title")] = string(R"(Element Swap Simplex move)");
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
	help_strings[string("mvEmpiricalTree")][string("title")] = string(R"(Move on an empirical tree distribution)");
	help_strings[string("mvEventTimeBeta")][string("name")] = string(R"(mvEventTimeBeta)");
	help_strings[string("mvEventTimeSlide")][string("name")] = string(R"(mvEventTimeSlide)");
	help_strings[string("mvFNPR")][string("name")] = string(R"(mvFNPR)");
	help_arrays[string("mvFossilTipTimeSlideUniform")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("mvFossilTipTimeSlideUniform")][string("description")] = string(R"(This moves either takes a specific fossil, or randomly picks a fossil, and then performs a sliding move on the tip age.)");
	help_strings[string("mvFossilTipTimeSlideUniform")][string("details")] = string(R"(This sliding move uses the possible minimum and maximum ages as reflection boundaries.
The maximum ages is computed either by its parents or the maximum age in the uncertainty of the fossil, which can be provided to the move or is taken from the taxon object.
The minimum ages is computed either by its oldest descendant (for sampled ancestors) or the minimum age in the uncertainty of the fossil, which can be provided to the move or is taken from the taxon object.)");
	help_strings[string("mvFossilTipTimeSlideUniform")][string("example")] = string(R"(
# Use a for loop to create a uniform distribution on the occurrence time for each fossil #
# The boundaries of the uniform distribution are specified in the tsv file #
fossils = fbd_tree.getFossils()
for(i in 1:fossils.size())
{
    t[i] := tmrca(fbd_tree, clade(fossils[i]))

    a[i] = fossils[i].getMinAge()
    b[i] = fossils[i].getMaxAge()

    F[i] ~ dnUniform(t[i] - b[i], t[i] - a[i])
    F[i].clamp( 0 )
    moves.append( mvFossilTipTimeUniform(fbd_tree, origin_time, min=a[i], max=b[i], tip=fossils[i], weight=5.0) )
    moves.append( mvFossilTipTimeSlideUniform(fbd_tree, origin_time, min=a[i], max=b[i], tip=fossils[i], weight=5.0) )
})");
	help_strings[string("mvFossilTipTimeSlideUniform")][string("name")] = string(R"(mvFossilTipTimeSlideUniform)");
	help_arrays[string("mvFossilTipTimeSlideUniform")][string("see_also")].push_back(string(R"(mvFossilTipTimeSlideUniform)"));
	help_strings[string("mvFossilTipTimeSlideUniform")][string("title")] = string(R"(Sliding move to change a fossil tip age)");
	help_arrays[string("mvFossilTipTimeUniform")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("mvFossilTipTimeUniform")][string("description")] = string(R"(This moves either takes a specific fossil, or randomly picks a fossil, and then draws the new ages randomly between the maximum and minimum ages.)");
	help_strings[string("mvFossilTipTimeUniform")][string("details")] = string(R"(The maximum ages is computed either by its parents or the maximum age in the uncertainty of the fossil, which can be provided to the move or is taken from the taxon object.
The minimum ages is computed either by its oldest descendant (for sampled ancestors) or the minimum age in the uncertainty of the fossil, which can be provided to the move or is taken from the taxon object.)");
	help_strings[string("mvFossilTipTimeUniform")][string("example")] = string(R"(
# Use a for loop to create a uniform distribution on the occurrence time for each fossil #
# The boundaries of the uniform distribution are specified in the tsv file #
fossils = fbd_tree.getFossils()
for(i in 1:fossils.size())
{
    t[i] := tmrca(fbd_tree, clade(fossils[i]))

    a[i] = fossils[i].getMinAge()
    b[i] = fossils[i].getMaxAge()

    F[i] ~ dnUniform(t[i] - b[i], t[i] - a[i])
    F[i].clamp( 0 )
    moves.append( mvFossilTipTimeUniform(fbd_tree, origin_time, min=a[i], max=b[i], tip=fossils[i], weight=5.0) )
    moves.append( mvFossilTipTimeSlideUniform(fbd_tree, origin_time, min=a[i], max=b[i], tip=fossils[i], weight=5.0) )
})");
	help_strings[string("mvFossilTipTimeUniform")][string("name")] = string(R"(mvFossilTipTimeUniform)");
	help_arrays[string("mvFossilTipTimeUniform")][string("see_also")].push_back(string(R"(mvFossilTipTimeSlideUniform)"));
	help_strings[string("mvFossilTipTimeUniform")][string("title")] = string(R"(Move to uniformly draw fossil tip ages)");
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
	help_arrays[string("mvIidPrior")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("mvIidPrior")][string("description")] = string(R"(This move proposes new values drawn from the prior.)");
	help_strings[string("mvIidPrior")][string("details")] = string(R"(Using this move, one actually gets an independence sampler as the proposal doesn't depend on the current state. The move calls redraw based on the distribution attached to the random variable.)");
	help_strings[string("mvIidPrior")][string("example")] = string(R"(x ~ dnUnif(0,10000)
moves[1] = mvIidPrior(x, weight=1.0)
monitors[1] = screenmonitor(printgen=1000, x)
mymodel = model(x)
mymcmc = mcmc(mymodel, monitors, moves)
mymcmc.run(generations=200000))");
	help_strings[string("mvIidPrior")][string("name")] = string(R"(mvIidPrio)");
	help_strings[string("mvIidPrior")][string("title")] = string(R"(Move to propose from prior)");
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
	help_arrays[string("mvRandomGeometricWalk")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("mvRandomGeometricWalk")][string("description")] = string(R"(A move that performs geometric random walk on an integer variable. The displacement of the random walk is drawn from a geometric distribution, mirrored for positive and negative steps.)");
	help_strings[string("mvRandomGeometricWalk")][string("example")] = string(R"(
p <- 0.8
x ~ dnGeom(p)

moves[1] = mvRandomGeometricWalk(x, weight=1.0)
monitors[1] = mvScreen(printgen=1000, x)

mymodel = model(p)
mymcmc = mcmc(mymodel, monitors, moves)
mymcmc.burnin(generations=20000,tuningInterval=100)
mymcmc.run(generations=200000))");
	help_strings[string("mvRandomGeometricWalk")][string("name")] = string(R"(mvRandomGeometricWalk)");
	help_arrays[string("mvRandomGeometricWalk")][string("see_also")].push_back(string(R"(mvRandomNaturalWalk)"));
	help_arrays[string("mvRandomGeometricWalk")][string("see_also")].push_back(string(R"(mvRandomIntegerWalk)"));
	help_strings[string("mvRandomGeometricWalk")][string("title")] = string(R"(Geometric random walk)");
	help_arrays[string("mvRandomIntegerWalk")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("mvRandomIntegerWalk")][string("description")] = string(R"(A move that performs random walk on an integer variable. The displacement of the random walk is exactly one step, either positive or negative.)");
	help_strings[string("mvRandomIntegerWalk")][string("example")] = string(R"(
p <- 0.8
x ~ dnGeom(p)

moves[1] = mvRandomIntegerWalk(x, weight=1.0)
monitors[1] = mvScreen(printgen=1000, x)

mymodel = model(p)
mymcmc = mcmc(mymodel, monitors, moves)
mymcmc.burnin(generations=20000,tuningInterval=100)
mymcmc.run(generations=200000))");
	help_strings[string("mvRandomIntegerWalk")][string("name")] = string(R"(
mvRandomIntegerWalk)");
	help_arrays[string("mvRandomIntegerWalk")][string("see_also")].push_back(string(R"(mvRandomNaturalWalk)"));
	help_arrays[string("mvRandomIntegerWalk")][string("see_also")].push_back(string(R"(mvRandomGeometricWalk)"));
	help_strings[string("mvRandomIntegerWalk")][string("title")] = string(R"(Random walk on integers)");
	help_arrays[string("mvRandomNaturalWalk")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("mvRandomNaturalWalk")][string("description")] = string(R"(A move that performs random walk on a natural number variable. The displacement of the random walk is exactly one step, either positive or negative.)");
	help_strings[string("mvRandomNaturalWalk")][string("example")] = string(R"(
p <- 0.8
x ~ dnGeom(p)

moves[1] = mvRandomNaturalWalk(x, weight=1.0)
monitors[1] = mvScreen(printgen=1000, x)

mymodel = model(p)
mymcmc = mcmc(mymodel, monitors, moves)
mymcmc.burnin(generations=20000,tuningInterval=100)
mymcmc.run(generations=200000))");
	help_strings[string("mvRandomNaturalWalk")][string("name")] = string(R"(
mvRandomNaturalWalk)");
	help_arrays[string("mvRandomNaturalWalk")][string("see_also")].push_back(string(R"(mvRandomIntegerWalk)"));
	help_arrays[string("mvRandomNaturalWalk")][string("see_also")].push_back(string(R"(mvRandomGeometricWalk)"));
	help_strings[string("mvRandomNaturalWalk")][string("title")] = string(R"(Random walk on natural numbers)");
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
	help_strings[string("mvSlice")][string("description")] = string(R"(Instead of using a fixed move size, `mvSlice` determines the size of a move proposal based on the current shape of the likelihood function.
This allows small moves to be proposed in certain parts of parameter space,
and large moves in other parts of the space, as appropriate.)");
	help_strings[string("mvSlice")][string("name")] = string(R"(Slice move)");
	help_arrays[string("mvSlice")][string("see_also")].push_back(string(R"(`mvSlide` and `mvScale` are possible alternatives where a fixed move size is desired.)"));
	help_strings[string("mvSlice")][string("title")] = string(R"(Propose a slice move)");
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
	help_strings[string("mvTreeScale")][string("name")] = string(R"(mvTreeScale)");
	help_strings[string("mvUPPAllocation")][string("name")] = string(R"(mvUPPAllocation)");
	help_strings[string("mvUpDownScale")][string("name")] = string(R"(mvUpDownScale)");
	help_strings[string("mvUpDownSlide")][string("name")] = string(R"(mvUpDownSlide)");
	help_strings[string("mvUpDownSlideBactrian")][string("name")] = string(R"(mvUpDownSlideBactrian)");
	help_strings[string("mvVectorBinarySwitch")][string("name")] = string(R"(mvVectorBinarySwitch)");
	help_arrays[string("mvVectorElementSwap")][string("authors")].push_back(string(R"(Sebastian Hoehna)"));
	help_strings[string("mvVectorElementSwap")][string("description")] = string(R"(Move that randomly picks a pair of elements in a vector on swaps the two with another.)");
	help_strings[string("mvVectorElementSwap")][string("name")] = string(R"(mvVectorElementSwap)");
	help_arrays[string("mvVectorElementSwap")][string("see_also")].push_back(string(R"(mvVectorBinarySwitch)"));
	help_strings[string("mvVectorElementSwap")][string("title")] = string(R"(Move to swap to elements in a vector)");
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
	help_strings[string("readVCF")][string("description")] = string(R"(Read VCF file into revbayes)");
	help_strings[string("readVCF")][string("details")] = string(R"(readVCF accepts two arguments to read in a VCF file. The first argument
specifies the relative or absolute file path to desired VCF file. The second
specifies type of data to be constructed (default binary). This function
only allows for 0, 1, and . characters in the VCF file.)");
	help_strings[string("readVCF")][string("example")] = string(R"(x <- readVCF("path/to/VCF/file", "DNA"))");
	help_strings[string("readVCF")][string("name")] = string(R"(readVCF)");
	help_strings[string("readVCF")][string("title")] = string(R"(Read VCF)");
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
	help_strings[string("setOption")][string("details")] = string(R"(Runtime options are used to personalize RevBayes and are stored on the local machine.
The currently available keys and their associated values are as follows:

    lineWidth=<integer>          Screen width when printing (in characters).
        
        DEFAULT: 160
    
    outputPrecision=<integer>    How many significant digits to print for the values of model graph nodes.
        
        DEFAULT: 7
    
    printNodeIndex=<true,false>  Print the node indices of a tree as annotations?
        
        DEFAULT: true
    
    useScaling=<true,false>      Should the partial likelihoods in continuous-time Markov chain (CTMC) models be scaled
                                 to avoid underflow?
        DEFAULT: true
    
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
        
        DEFAULT: 0)");
	help_strings[string("setOption")][string("example")] = string(R"(# compute the absolute value of a real number
getOption("linewidth")

# let us set the linewidth to a new value
setOption("linewidth", 200)

# now let's check what the value is
getOption("linewidth"))");
	help_strings[string("setOption")][string("name")] = string(R"(setOption)");
	help_arrays[string("setOption")][string("see_also")].push_back(string(R"(getOption)"));
	help_strings[string("setOption")][string("title")] = string(R"(Set a global RevBayes option)");
	help_strings[string("setValue")][string("description")] = string(R"(`x.setValue(value)` sets the value of the stochastic variable `x` to `value`.)");
	help_strings[string("setValue")][string("details")] = string(R"(`x.setValue()` allows calculations to be evaluated at a given value of `x`, whilst allowing the value
of `x` to be redrawn, or to vary during MCMC.

`.setValue()` allows an MCMC run to be initialized with plausible values,
which can expedite convergence when priors are broad
([example](https://revbayes.github.io/tutorials/divrate/branch_specific.html#specifying-the-model)),
and can be useful when [debugging MCMC runs](https://revbayes.github.io/tutorials/mcmc_troubleshooting/#starting-values).)");
	help_strings[string("setValue")][string("example")] = string(R"(x ~ dnNormal(1, 1)

# Evaluate P(x) at x = 1
x.setValue(1)
x.probability()

# Modify the observed value of x
x.redraw()
x.probability()

# Initialize an MCMC run with a specific value
x.setValue(40000)
mcmc(model(x), [mnScreen(x)], [mvScale(x)]).run(generations = 5))");
	help_strings[string("setValue")][string("name")] = string(R"(Set value)");
	help_arrays[string("setValue")][string("see_also")].push_back(string(R"(clamp)"));
	help_strings[string("setValue")][string("title")] = string(R"(Set the value of a stochastic variable)");
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
	help_strings[string("srGelmanRubin")][string("description")] = string(R"(Allow an MCMC run to terminate once the specified criterion has been met.
The Gelman–Rubin rule compares the variance between runs with the variance within runs; its value tends to unity (1) as runs converge.  It is widely referred to as the "potential scale reduction factor" (PSRF).)");
	help_strings[string("srGelmanRubin")][string("details")] = string(R"(Because the statistic is defined by comparing the variation between different runs to the variance within each run, it can only be calculated when multiple independent runs are performed, by setting the `nruns` argument to `mcmc` or `mcmcmc` to a value greater than one.)");
	help_strings[string("srGelmanRubin")][string("example")] = string(R"(```
# Binomial example: estimate success probability given 7 successes out of 20 trials
r ~ dnExp(10)
p := Probability(ifelse(r < 1, r, 1))
n <- 20
k ~ dnBinomial(n, p)
k.clamp(7)
mymodel = model(k)

moves = VectorMoves()
moves.append( mvSlide(r, delta=0.1, weight=1) )

paramFile = "parameters.log"

monitors = VectorMonitors()
monitors.append( mnModel(filename=paramFile, printgen=100, p) )

# Stop when the potential scale reduction factor falls below 1.01
stopping_rules[1] = srGelmanRubin(1.01, file = paramFile, freq = 1000)

# Create the MCMC object
mymcmc = mcmc(mymodel, monitors, moves, nruns = 2)

# Begin the MCMC run
mymcmc.run(rules = stopping_rules)
```)");
	help_strings[string("srGelmanRubin")][string("name")] = string(R"(srGelmanRubin)");
	help_references[string("srGelmanRubin")].push_back(RbHelpReference(R"(Gelman, A; Rubin, D.B. (1992). Inference from Iterative Simulation Using Multiple Sequences. Statistical Science. 7 (4): 457–472)",R"(10.1214/ss/1177011136)",R"()"));
	help_references[string("srGelmanRubin")].push_back(RbHelpReference(R"(Vats, D.; Knudson, C. Revisiting the Gelman–Rubin Diagnostic. Statist. Sci. 36 (4) 518 )",R"()",R"()"));
	help_references[string("srGelmanRubin")].push_back(RbHelpReference(R"()",R"(10.1214/20-STS812 )",R"()"));
	help_arrays[string("srGelmanRubin")][string("see_also")].push_back(string(R"()"));
	help_arrays[string("srGelmanRubin")][string("see_also")].push_back(string(R"(- Tutorial on [convergence assessment](https://revbayes.github.io/tutorials/convergence/))"));
	help_strings[string("srGelmanRubin")][string("title")] = string(R"(Gelman–Rubin (PSRF) stopping rule)");
	help_arrays[string("srGeweke")][string("authors")].push_back(string(R"(Incorporates text by Martyn Plummer)"));
	help_strings[string("srGeweke")][string("description")] = string(R"(Allow an MCMC run to terminate once the specified criterion has been met.

Geweke (1992) proposed a convergence diagnostic for Markov chains based on a test for equality of the means of the first and last part of a Markov chain (by default the first 10% and the last 50%). If the samples are drawn from the stationary distribution of the chain, the two means are equal and Geweke's statistic has an asymptotically standard normal distribution.

The test statistic is a standard Z-score: the difference between the two sample means divided by its estimated standard error. The standard error is estimated from the spectral density at zero and so takes into account any autocorrelation.

The Z-score is calculated under the assumption that the two parts of the chain are asymptotically independent, which requires that the sum of `frac1` and `frac2` be strictly less than 1.)");
	help_strings[string("srGeweke")][string("example")] = string(R"(```
# Binomial example: estimate success probability given 7 successes out of 20 trials
r ~ dnExp(10)
p := Probability(ifelse(r < 1, r, 1))
n <- 20
k ~ dnBinomial(n, p)
k.clamp(7)
mymodel = model(k)

moves = VectorMoves()
moves.append( mvSlide(r, delta=0.1, weight=1) )

paramFile = "parameters.log"

monitors = VectorMonitors()
monitors.append( mnModel(filename=paramFile, printgen=100, p) )

# Stop when the Geweke test statistic becomes significant at alpha = 0.001
stopping_rules[1] = srGeweke( prob=0.001, file=paramFile, freq=10000 )

# Create the MCMC object
mymcmc = mcmc( mymodel, monitors, moves )

# Begin the MCMC run
mymcmc.run( rules = stopping_rules )
```)");
	help_strings[string("srGeweke")][string("name")] = string(R"(srGeweke)");
	help_references[string("srGeweke")].push_back(RbHelpReference(R"(Geweke, J. Evaluating the accuracy of sampling-based approaches to calculating posterior moments.  In Bayesian Statistics 4 (ed JM Bernado, JO Berger, AP Dawid and AFM Smith).   Clarendon Press, Oxford, UK. )",R"()",R"()"));
	help_arrays[string("srGeweke")][string("see_also")].push_back(string(R"()"));
	help_arrays[string("srGeweke")][string("see_also")].push_back(string(R"(- Tutorial on [convergence assessment](https://revbayes.github.io/tutorials/convergence/))"));
	help_strings[string("srGeweke")][string("title")] = string(R"(Geweke stopping rule)");
	help_strings[string("srMaxIteration")][string("description")] = string(R"(Cause an MCMC run to terminate once the specified number of iterations have been performed.
This function would typically be used alongside other stopping rules)");
	help_strings[string("srMaxIteration")][string("example")] = string(R"(```
# Binomial example: estimate success probability given 7 successes out of 20 trials
r ~ dnExp(10)
p := Probability(ifelse(r < 1, r, 1))
n <- 20
k ~ dnBinomial(n, p)
k.clamp(7)
mymodel = model(k)

moves = VectorMoves()
moves.append( mvSlide(r, delta=0.1, weight=1) )

paramFile = "parameters.log"

monitors = VectorMonitors()
monitors.append( mnModel(filename=paramFile, printgen=100, p) )

# Stop when 1000 iterations have been completed
stopping_rules[1] = srMaxIteration(1000)

# Create the MCMC object
mymcmc = mcmc(mymodel, monitors, moves)

# Begin the MCMC run
mymcmc.run(rules = stopping_rules)
```)");
	help_strings[string("srMaxIteration")][string("name")] = string(R"(srMaxIteration)");
	help_strings[string("srMaxIteration")][string("title")] = string(R"(Maximum iteration stopping rule)");
	help_strings[string("srMaxTime")][string("description")] = string(R"(Cause an MCMC run to terminate once the specified time has elapsed.)");
	help_strings[string("srMaxTime")][string("example")] = string(R"(```

# Binomial example: estimate success probability given 7 successes out of 20 trials
r ~ dnExp(10)
p := Probability(ifelse(r < 1, r, 1))
n <- 20
k ~ dnBinomial(n, p)
k.clamp(7)
mymodel = model(k)

moves = VectorMoves()
moves.append( mvSlide(r, delta=0.1, weight=1) )

paramFile = "parameters.log"

monitors = VectorMonitors()
monitors.append( mnModel(filename=paramFile, printgen=100, p) )

# Stop when the five seconds have elapsed
stopping_rules[1] = srMaxTime(5, "seconds")

# Create the MCMC object
mymcmc = mcmc(mymodel, monitors, moves)

# Begin the MCMC run
mymcmc.run(rules = stopping_rules)
```)");
	help_strings[string("srMaxTime")][string("name")] = string(R"(srMaxTime)");
	help_strings[string("srMaxTime")][string("title")] = string(R"(Maximum time stopping rule)");
	help_arrays[string("srMinESS")][string("authors")].push_back(string(R"(ESS explanation adapted from Luiza Fabreti and Sebastian Höhna's [tutorial](https://revbayes.github.io/tutorials/convergence/))"));
	help_strings[string("srMinESS")][string("description")] = string(R"(Allow an MCMC run to terminate once the specified criterion has been met.)");
	help_strings[string("srMinESS")][string("details")] = string(R"(The Effective Sample Size (ESS) is the number of independent samples generated by a MCMC sampler.
The ESS takes into account the correlation between samples within a chain.
Low ESS values represent high autocorrelation in the chain.
If the autocorrelation is higher, then the uncertainty in our estimates is also higher.

The MCMC run will terminate once all parameters in every log file meet the ESS
threshold.  As such, performing additional runs will not decrease the number
of generations required to meet the ESS threshold – even though it will increase
the number of indepedent samples in the final, pooled posterior sample.)");
	help_strings[string("srMinESS")][string("example")] = string(R"(```
# Binomial example: estimate success probability given 7 successes out of 20 trials
r ~ dnExp(10)
p := Probability(ifelse(r < 1, r, 1))
n <- 20
k ~ dnBinomial(n, p)
k.clamp(7)
mymodel = model(k)

moves = VectorMoves()
moves.append( mvSlide(r, delta=0.1, weight=1) )

paramFile = "parameters.log"

monitors = VectorMonitors()
monitors.append( mnModel(filename=paramFile, printgen=100, p) )

# Stop when all monitored parameters have attained an estimated sample size of 50
stopping_rules[1] = srMinESS(50, file = paramFile, freq = 1000)

# Create the MCMC object
mymcmc = mcmc(mymodel, monitors, moves)

# Begin the MCMC run
mymcmc.run(rules = stopping_rules)
```)");
	help_strings[string("srMinESS")][string("name")] = string(R"(srMinESS)");
	help_arrays[string("srMinESS")][string("see_also")].push_back(string(R"()"));
	help_arrays[string("srMinESS")][string("see_also")].push_back(string(R"(- The tutorial on [convergence assessment](https://revbayes.github.io/tutorials/convergence/) contains a discusson on the calculation and interpretation of the ESS diagnostic.)"));
	help_strings[string("srMinESS")][string("title")] = string(R"(Estimated sample size stopping rule)");
	help_strings[string("srStationarity")][string("description")] = string(R"(Allow an MCMC run to terminate once the specified criterion has been met.
An MCMC sample can be considered stationary once its mean, variance and autocorrelation structure do not change over time.)");
	help_strings[string("srStationarity")][string("details")] = string(R"(Because the statistic is defined by comparing different runs, it can only be calculated when multiple independent runs are performed, by setting the `nruns` argument to `mcmc` or `mcmcmc` to a value greater than one.)");
	help_strings[string("srStationarity")][string("example")] = string(R"(```
# Binomial example: estimate success probability given 7 successes out of 20 trials
r ~ dnExp(10)
p := Probability(ifelse(r < 1, r, 1))
n <- 20
k ~ dnBinomial(n, p)
k.clamp(7)
mymodel = model(k)

moves = VectorMoves()
moves.append( mvSlide(r, delta=0.1, weight=1) )

paramFile = "parameters.log"

monitors = VectorMonitors()
monitors.append( mnModel(filename=paramFile, printgen=100, p) )

# Stop when stationarity has been attained at confidence level gamma = 0.25
stopping_rules[1] = srStationarity(prob = 0.25, file = paramFile, freq = 1000)

# Create the MCMC object
mymcmc = mcmc(mymodel, monitors, moves, nruns = 2)

# Begin the MCMC run
mymcmc.run(rules = stopping_rules)
```)");
	help_strings[string("srStationarity")][string("name")] = string(R"(srStationarity)");
	help_references[string("srStationarity")].push_back(RbHelpReference(R"(Hill, S.D. and Spall, J.C. 2011. Stationarity and Convergence of the Metropolis-Hastings Algorithm: Insights into Theoretical Aspects. IEEE Control Systems Magazine 39.)",R"(10.1109/MCS.2018.2876959)",R"()"));
	help_arrays[string("srStationarity")][string("see_also")].push_back(string(R"(- Tutorial on [convergence assessment](https://revbayes.github.io/tutorials/convergence/))"));
	help_strings[string("srStationarity")][string("title")] = string(R"(Stationarity stopping rule)");
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
	help_strings[string("sum")][string("description")] = string(R"(Sums all members of a vector of type `Real`, `RealPos`, `Integer`, or `Natural`)");
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
	help_strings[string("time")][string("description")] = string(R"(Get the current system time.)");
	help_strings[string("time")][string("details")] = string(R"(
"year" reports the current year (e.g. 2000).

"day" returns the index of the day in the year (e.g. Jan 1 = 1; Feb 1 = 32).

"(milli)seconds" returns the number of (milli)seconds that have elapsed since midnight.

"fromBeginning", the default, returns the number of milliseconds that have elapsed since 1400-Jan-01 00:00:00, the earliest representable date in the boost library's implementation of the Gregorian date system.)");
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
	help_arrays[string("tnExp")][string("authors")].push_back(string(R"(Ben Redelings)"));
	help_strings[string("tnExp")][string("description")] = string(R"(Exp-transforms a given distribution.)");
	help_strings[string("tnExp")][string("details")] = string(R"(If X ~ dist then tnExp(dist) is the distribution of exp(X).

The distribution `dist` can be either univariate (dnNormal) or
multivariate (dnMultivariateNormal).

This turns out to be the same as dnLog(dist), which provides a distribution
that has distribution `dist` on the log-scale.)");
	help_strings[string("tnExp")][string("example")] = string(R"(x ~ tnExp(dnNormal(0,1))          # Draw from the log-Normal distribution
x ~ dnNormal(0,1) |> tnExp()      # Expressed using pipes.
x ~ dnLognormal(0,1)              # This is equivalent.
y ~ dnNormal(0,1)
x := exp(y)                       # This is also equivalent.

x ~ tnExp(dnGamma(2,3))           # There is no equivalent for this.
x ~ dnIID(10,tnExp(dnGamma(2,3))) # Draw 10 log-Gamma(2,3) random variables.

mu = [1.0, 2.0, 3.0, 4.0]
Sigma ~ dnWishart(df=4, kappa=2, dim=4)
x ~ dnMultivariateNormal(mu,Sigma) |> tnExp())");
	help_strings[string("tnExp")][string("name")] = string(R"(tnExp)");
	help_arrays[string("tnExp")][string("see_also")].push_back(string(R"(tnLog, tnLogit, tnInvlogit)"));
	help_strings[string("tnExp")][string("title")] = string(R"(Exp-transformed distribution)");
	help_arrays[string("tnInvlogit")][string("authors")].push_back(string(R"(Ben Redelings)"));
	help_strings[string("tnInvlogit")][string("description")] = string(R"(Invlogit-transforms a given distribution.)");
	help_strings[string("tnInvlogit")][string("details")] = string(R"(If X ~ dist then tnInvlogit(dist) is the distribution of exp(X)/(1+exp(X)).
The inverse logit function is also called the logistic function.

The distribution `dist` can be either univariate (dnNormal) or
multivariate (MultivariateNormal).)");
	help_strings[string("tnInvlogit")][string("example")] = string(R"(p ~ tnInvlogit(dnNormal(0,1))      # The inverse-logit of a Normal random variable.
p ~ dnNormal(0,1) |> tnInvlogit()  # Expressed using pipes.

x ~ dnNormal(0,1)
p := invlogit(x)                   # Expressed as a deterministic function of the log-odds.

ps ~ dnIID(4,dnNormal(0,1)) |> tnInvlogit()

mu = [1.0, 2.0, 3.0, 4.0]
Sigma ~ dnWishart(df=4, kappa=2, dim=4)
x ~ dnMultivariateNormal(mu,Sigma) |> tnInvlogit())");
	help_strings[string("tnInvlogit")][string("name")] = string(R"(tnInvlogit)");
	help_arrays[string("tnInvlogit")][string("see_also")].push_back(string(R"(logistic, tnExp, tnLog, tnLogit)"));
	help_strings[string("tnInvlogit")][string("title")] = string(R"(Invlogit-transformed distribution)");
	help_arrays[string("tnLog")][string("authors")].push_back(string(R"(Ben Redelings)"));
	help_strings[string("tnLog")][string("description")] = string(R"(Log-transforms a given distribution.)");
	help_strings[string("tnLog")][string("details")] = string(R"(If X ~ dist then tnLog(dist) is the distribution of log(X).

The distribution `dist` can be either univariate (dnExponential) or
multivariate (dnDirichlet).

This is NOT the same as dnLog(dist), which provides a distribution
that has distribution `dist` on the log-scale.)");
	help_strings[string("tnLog")][string("example")] = string(R"(x ~ tnLog(dnExponential(1))       # The log of an Exponential random variable.
x ~ dnExponential(1) |> tnLog()   # Expressed using pipes.

y ~ dnExponential(1)
x := log(y)                       # This is also equivalent.

x ~ dnDirichlet([1,1,1,1]) |> tnLog())");
	help_strings[string("tnLog")][string("name")] = string(R"(tnLog)");
	help_arrays[string("tnLog")][string("see_also")].push_back(string(R"(tnExp, tnLogit, tnInvlogit)"));
	help_strings[string("tnLog")][string("title")] = string(R"(Log-transformed distribution)");
	help_arrays[string("tnLogit")][string("authors")].push_back(string(R"(Ben Redelings)"));
	help_strings[string("tnLogit")][string("description")] = string(R"(Logit-transforms a given distribution.)");
	help_strings[string("tnLogit")][string("details")] = string(R"(If P ~ dist then tnLogit(dist) is the distribution of log(P/(1-P)).

The distribution `dist` can be either univariate (dnBeta) or
multivariate (dnDirichlet).)");
	help_strings[string("tnLogit")][string("example")] = string(R"(x ~ tnLogit(dnBeta(1,2))         # The log-odds of an Beta random variable.
x ~ dnBeta(1,2)|> tnLogit()      # Expressed using pipes.

p ~ dnBeta(1,2)
x := logit(p)                    # Expressed as a deterministic function of the probability.

xs ~ dnDirichlet([1,1,1,1]) |> tnLogit())");
	help_strings[string("tnLogit")][string("name")] = string(R"(tnLogit)");
	help_arrays[string("tnLogit")][string("see_also")].push_back(string(R"(logit, tnExp, tnLog, tnInvlogitit)"));
	help_strings[string("tnLogit")][string("title")] = string(R"(Logit-transformed distribution)");
	help_arrays[string("tnScale")][string("authors")].push_back(string(R"(Ben Redelings)"));
	help_strings[string("tnScale")][string("description")] = string(R"(Rescales a given distribution.)");
	help_strings[string("tnScale")][string("details")] = string(R"(If X ~ dist then tnScale(dist, lambda) is the distribution of X * lambda)");
	help_strings[string("tnScale")][string("example")] = string(R"(x ~ tnScale(dExponential(1),2)       # An Exponential(rate=0.5) random variable.
x ~ dnExponential(1) |> tnScale(2)   # Expressed using pipes.)");
	help_strings[string("tnScale")][string("name")] = string(R"(tnScale)");
	help_arrays[string("tnScale")][string("see_also")].push_back(string(R"(tnScale)"));
	help_strings[string("tnScale")][string("title")] = string(R"(A scaled distribution)");
	help_arrays[string("tnShift")][string("authors")].push_back(string(R"(Ben Redelings)"));
	help_strings[string("tnShift")][string("description")] = string(R"(Shifts a given distribution.)");
	help_strings[string("tnShift")][string("details")] = string(R"(If X ~ dist then tnShift(dist, d) is the distribution of X + d)");
	help_strings[string("tnShift")][string("example")] = string(R"(x ~ tnShift(dExponential(1),2)       # An exponential variable starting at 2.
x ~ dnExponential(1) |> tnShift(2)   # Expressed using pipes.)");
	help_strings[string("tnShift")][string("name")] = string(R"(tnShift)");
	help_arrays[string("tnShift")][string("see_also")].push_back(string(R"(tnScale)"));
	help_strings[string("tnShift")][string("title")] = string(R"(A shifted distribution)");
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
	help_arrays[string("vectorFlatten")][string("authors")].push_back(string(R"(Michael Landis)"));
	help_strings[string("vectorFlatten")][string("description")] = string(R"(Flatten a vector to one dimension.)");
	help_strings[string("vectorFlatten")][string("details")] = string(R"(This function accepts a two-dimensional vector as an argument and flattens it
to one dimension.)");
	help_strings[string("vectorFlatten")][string("example")] = string(R"(# Define Vector
x <- v([1, 2], [3, 4])
# Flatten
x_flat <- vectorFlatten(x)
x_flat
[1, 2, 3, 4])");
	help_strings[string("vectorFlatten")][string("name")] = string(R"(vectorFlatten)");
	help_arrays[string("vectorFlatten")][string("see_also")].push_back(string(R"(v)"));
	help_arrays[string("vectorFlatten")][string("see_also")].push_back(string(R"(RealPos)"));
	help_strings[string("vectorFlatten")][string("title")] = string(R"(Vector Flatten)");
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
