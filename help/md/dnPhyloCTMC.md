## name
`dnPhyloCTMC`: Distribution of a phylogenetic continuous-time Markov chain

## title
The parameters of a phylogenetic model – a tree topology with branch lengths, a substitution model that describes how observations evolve over the tree, etc. – collectively form a distribution called the _phylogenetic continuous-time Markov chain_.

## description
dnPhyloCTMC gives the probability distribution of tip labels on a phylogenetic tree given an phylogenetic continuous-time Markov chain model.

## details

The likelihood of observed tip labels (specified via a clamped `AbstractHomologousDiscreteCharacterData` object) is computed using Felsenstein's pruning algorithm, with partial likelihoods stored for each branch of the tree. It is automatically outputted in the `Likelihood` column of the `mnFile()` and `mnScreen()` monitors (which can be suppressed with `likelihood = FALSE`).

## authors
## see_also
- Tutorial on [graphical models](https://revbayes.github.io/tutorials/intro/graph_models)

- Tutorial on [specifying a phylogenetic continuous-time Markov chain](https://revbayes.github.io/tutorials/ctmc/) model

## example

```rb
# Create stochastic node with the tip distribution given by the topology `tree` with the 
# rate matrix parameter q_matrix
x ~ dnPhyloCTMC(tree = tree, Q = q_matrix)

# Clamp observed characters to the node
chars <- readDiscreteCharacterData("myData.nex")
x.clamp(chars)

# Calculate the probability of the observed characters under the given distribution
x.lnProbability()
```

## Methods

MixtureLikelihoods = MixtureLikelihoods ()

[] = [] (Natural<any> index)

addMissingTaxa = addMissingTaxa (String|Taxon|String[]|Taxon[]<any>...

applyMissingSitesMask = applyMissingSitesMask (AbstractHomologousDi...

chartype = chartype ()

clamp = clamp (AbstractHomologousDiscreteCharacterData)[https://revbayes.github.io/documentation/AbstractHomologousDiscreteCharacterData.html])<any> x): Clamp the draw from this distribution to a `AbstractHomologousDiscreteCharacterData` object.

computeMultinomialProfileLikelihood = computeMultinomialProfileLike...

computeSiteFrequencySpectrum = computeSiteFrequencySpectrum (Bool<a...

computeStateFrequencies = computeStateFrequencies ()

excludeAll = excludeAll ()

excludeCharacter = excludeCharacter (Natural<any> pos)

excludeCharacter = excludeCharacter (Natural[]<any> )

excludeMissingSites = excludeMissingSites ()

excludeTaxa = excludeTaxa (String|Taxon<any> taxon)

excludeTaxa = excludeTaxa (String[]|Taxon[]<any> taxa)

expandCharacters = expandCharacters (Natural<any> factor)

filename = filename ()

getEmpiricalBaseFrequencies = getEmpiricalBaseFrequencies ()

getIncludedCharacterIndices = getIncludedCharacterIndices ()

getInvariantSiteIndices = getInvariantSiteIndices (Bool<any> exclud...

getNumInvariantSites = getNumInvariantSites (Bool<any> excludeAmbig...

getNumStatesVector = getNumStatesVector ()

getPairwiseDifference = getPairwiseDifference (Bool<any> excludeAmb...

getStateDescriptions = getStateDescriptions ()

includeAll = includeAll ()

includeCharacter = includeCharacter (Natural<any> )

includeCharacter = includeCharacter (Natural[]<any> )

includeTaxa = includeTaxa (String<any> name)

includeTaxa = includeTaxa (String[]<any> names)

integrateOut = integrateOut (Bool<any> x)

isHomologous = isHomologous ()

isResolved = isResolved (Natural<any> taxonIndex,...

isSequenceMissing = isSequenceMissing (String<any> name)

lnMixtureLikelihoods = lnMixtureLikelihoods ()

lnProbability = lnProbability ()

maxGcContent = maxGcContent (Bool<any> excludeAmbiguous)

maxInvariableBlockLength = maxInvariableBlockLength (Bool<any> excl...

maxPairwiseDifference = maxPairwiseDifference (Bool<any> excludeAmb...

maxStates = maxStates ()

maxVariableBlockLength = maxVariableBlockLength (Bool<any> excludeA...

meanGcContent = meanGcContent (Bool<any> excludeAmbiguous)

meanGcContentByCodonPosition = meanGcContentByCodonPosition (Natura...

methods = methods (): List all available methods

minGcContent = minGcContent (Bool<any> excludeAmbiguous)

minPairwiseDifference = minPairwiseDifference (Bool<any> excludeAmb...

names = names ()

nchar = nchar ()

ntaxa = ntaxa ()

numInvariableBlocks = numInvariableBlocks (Bool<any> excludeAmbiguo...

numTaxaMissingSequence = numTaxaMissingSequence (Probability<any> x)

percentageMissing = percentageMissing (String<any> name)

probability = probability ()

redraw = redraw ()

removeExcludedCharacters = removeExcludedCharacters ()

removeTaxa = removeTaxa (String<any> name)

removeTaxa = removeTaxa (String[]<any> names)

replaceRandomSitesByMissingData = replaceRandomSitesByMissingData (...

setCodonPartition = setCodonPartition (Natural<any> )

setCodonPartition = setCodonPartition (Natural[]<any> )

setHomeologPhase = setHomeologPhase (String<any> data_name,...

setNumStatesPartition = setNumStatesPartition (Natural<any> )

setTaxonName = setTaxonName (String<any> current,...

setTaxonObject = setTaxonObject (String<any> current,...

setValue = setValue (AbstractHomologousDiscreteCharacterData<any> x)

show = show ()

siteLikelihoods = siteLikelihoods ()

siteMixtureLikelihoods = siteMixtureLikelihoods ()

siteRateLikelihoods = siteRateLikelihoods ()

siteRates = siteRates (String<any> estimateMethod {valid options: "...

size = size ()

taxa = taxa ()

taxonIndex = taxonIndex (String<any> name)

translateCharacters = translateCharacters (String<any> type)

unclamp = unclamp ()

varGcContent = varGcContent (Bool<any> excludeAmbiguous)

varGcContentByCodonPosition = varGcContentByCodonPosition (Natural<...


## references
- citation: Felsenstein J., 1973. Maximum Likelihood and Minimum-Steps Methods for Estimating Evolutionary Trees from Data on Discrete Characters. Systematic Biology 22:3, 240--249
  doi = 10.1093/sysbio/22.3.240
- citation: Felsenstein, J. (1981). Evolutionary trees from DNA sequences: A maximum likelihood approach. Journal of Molecular Evolution. 17 (6): 368–376.
  doi : 10.1007/BF01734359
- citation: Hohna, S., Landis, M.J. and Heath, T.A. 2017. Phylogenetic inference using `RevBayes`. Curr. Protoc. Bioinform.
57:6.16.1-6.16.34.
  doi: 10.1002/cpbi.22
  url: null


