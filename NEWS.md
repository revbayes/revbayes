# RevBayes 1.4.0 (unreleased)

## Backwards-incompatible changes

## Features

## Bug fixes

## Documentation improvements

## Infrastructure

# RevBayes 1.3.1 (Aug 7, 2025)

## Backwards-incompatible changes
  * Change the distribution name of the sampled-speciation birth-death model from `dnSBBDP` to `dnSSBDP` (#714).

## Features
  * Methods / arguments
      - Add an `.nbranches( )` method to `Tree` (#759).
      - Add a method for calculating the standard error of marginal likelihood estimates (#779).
      - Add a nanoseconds option to `time( )` (#796).
  * Distributions
      - Add a topologically constrained distribution for the state-dependent speciation and extinction process (#693).
      - Allow nonzero `rho` values in FBD analyses without extant tips (#788).

## Bug fixes
  * Member procedures
      - Fix `.getTopologyFrequency( )` behavior for tree traces (#750).
      - Fix `.resolveMultifurcations( )` behavior in edge cases (#771, #780).
  * Nexus I/O
      - Make sure `writeNexus()` can handle `BranchLengthTree` vectors (#758).
      - Fix parsing of the `STATELABELS` block in Nexus files (#760).
      - Allow reading Nexus files with an `ASSUMPTIONS` block (#764).
  * MPI / MC^3
      - Fix segfault in MPI + MC^3 (#763).
      - Fix logging and progress monitoring for `.runOneStone()` + MPI (#782).
  * Moves
      - Fix argument name mismatch in `mvUpDownSlide` (#766).
      - Fix `mvEllipticalSliceSamplingSimple` to not require strict equality of likelihoods (#773).
      - Reject tree moves if there are no nodes for them to operate on (#798, #802).
  * Misc
      - Fix and extend TensorPhylo, stochastic mapping, and biogeography models (#693).
      - Synchronize dirty flagging for nodes and transition probability matrices (#773).
      - Fix caching to avoid crashes when using `dnReversibleJumpMixture` (#773).
      - Fix posterior predictive simulations (#778).
  * Partial
      - Avoid subscripting by `-1` in `dnAutocorrelatedEvent` (#773).

## Documentation improvements
  * Basic data types: `Integer`, `Natural`, `Probability`, `Real` (#714).
  * Birth-death models: `dnCBDSP`, `dnCDBDP`, `dnSSBDP` (#714).
  * Marginal likelihood estimators: `pathSampler`, `steppingStoneSampler` (#779).

## Infrastructure
  * Update DLLs to allow running RevStudio on Windows (#751).
  * Don't print operator summaries in coalescent tests, and only run 1 generation (#770).
  * Enable testing with assertions to detect out-of-bounds indices in `std::vector` (#773).
  * Don't skip remaining scripts in a test if one fails, and combine error message from multiple failing scripts in a test (#773).
  * Allow specifying a required exit code to allow/require failure (#773).
  * Fix check for missing output files (#773).
  * Make test runner output clearer and more informative (#773).

## New contributors
  * @prilau made their first contribution in #796.

# RevBayes 1.3.0 (May 2, 2025)

## Backwards-incompatible changes
  * Change the name of the tuning argument of `mvUpDownSlide` from `lambda` to `delta`.
  * Change the name of the first argument of `Tree.reroot( )` from `clade` to `outgroup`.
  * In `powerPosterior( )`, specifying `cats=N` sets up N analyses numbered 1--N, rather than N+1 analyses numbered 0--N.

## Features
  * Interface
      - Print convergence statistics during MCMC runs with stopping rules if `verbose=2` (#645, #657, #739).
      - Better progress monitoring for power posterior analyses (#673).
  * Methods / arguments
      - Add member procedure for randomly resolving multifurcations (#641).
      - Add an `nruns` argument to `readTrace( )` (#736).
  * Functions
      - Add function for minimum branch length time-scaling (#641).
      - Add trigonometric functions (#647, #648).

## Bug fixes
  * Checkpointing
      - Fix segfault when attempting to resume an MPI analysis (#646).
      - Fix mismatch of chain states and chain heats when resuming an MC^3 analysis (#670).
      - Don't skip iterations when logging a resumed MC^3 analysis (#727).
  * MPI / MC^3
      - Fix NaN heats in MC^3 (#661).
      - Make sure MPI runs the specified number of independent replicates (#662).
  * Model graph
      - Only convert deterministic nodes to constant nodes if they really are constant (#678).
      - Fix crash in user functions (#687).
      - Make sure indexing by a deterministic variable creates a deterministic node (#719).
      - Only type-convert variables based on value if they are constant (#721).
  * Moves
      - Fix `mvRateAgeBetaShift` (#688).
      - Allow moves to propose out-of-range indices (#709).
      - Allow moves with low probability ratio and high Jacobian/Hastings ratio (#728).
  * Member procedures
      - Fix `.dropTip( )` behavior for multifurcations (#640).
      - Fix `.getStateDescriptions( )` behavior when subsetting by state space (#696).
      - Fix `.reroot( )` behavior when specifying an outgroup by its name string (#742).
  * Misc
      - Fix stochastic character mapping when there are excluded characters (#636).
      - Fix handling of large integers on Windows (#708).
      - Don't ignore negative clade constraints (#711).
      - Allow for machine uncertainty in `dnEpisodicBirthDeath` with empirical sampling (#713).
      - Fix partial likelihood caching in `dnPhyloCTMC` (#729).
      - Allow reading multiple trees from a string (`readTrees(text = ...)`) (#735).
  * Partial
      - Prevent some instances of clamped values from being modified (#600).

## Documentation improvements
  * `Simplex` (#606).
  * Moves
      - `mvBetaSimplex`, `mvDirichletSimplex`, `mvElementSwapSimplex` (#606).
      - `mvDPPValueBetaSimplex`, `mvDPPValueScaling`, `mvDPPValueSliding` (#666).
      - `mvNNI`, `mvSPR`, `mvScale`, `mvScaleBactrian`, `mvSlide`, `mvSlideBactrian`, `mvTreeScale`, `mvUpDownScale`, `mvUpDownSlide` (#683).
  * `sin` (#648, #683).
  * Substitution models
      - `fnJC` (#649).
      - `fnF81`, `fnGTR`, `fnHKY`, `fnK80`, `fnK81`, `fnT92`, `fnTrN` (#653, #656).
      - `fnLG`, `fnWAG` (#655).
      - `fnFreeK` (#667).
      - `fnCovarion` (#704).
      - `fnFreeBinary` (#706).
  * `exp` (#653).
  * Complete and standardize documentation for `dnBernoulli`, `dnBeta`, `dnBimodalLognormal`, `dnBimodalNormal`, `dnBinomial`, `dnCategorical`, `dnCauchy`, `dnChisq`, `dnDirichlet`, `dnExponential`, `dnGamma` (#663).
  * `dnUniformInteger`, `fnReadVCF`, `vectorFlatten` (#667).
  * Expand documentation for `powerPosterior` (#673).
  * `matrix`, `var` (#681).
  * `floor`, `mnModel`, `reverse`, `sinh` (#683).
  * `power`, `write` (#690).
  * Analysis output types and I/O functions
      - `readTrace` (#702, #723, #736).
      - `readTreeTrace` (#703, #722, #723, #736).
      - `Trace`, `TraceTree` (#723).
  * Discretization functions: `fnDiscretizeDistribution`, `fnDiscretizeGamma` (#707).

## Infrastructure
  * Update validation tests (#473).
  * Fix continuous-integration builds by dropping openlibm (#644).
  * Change build files to streamline updates to the help database (#659).
  * Make tutorial tests more flexible (#674).
  * Make sure the test runner can handle tutorial checkpoint tests (#676).
  * Add integration tests for revised coalescent classes (#689).
  * Change website submodule to pull from source rather than master (#692).

## New contributors
  * @sigibrock made their first contribution in #648.
  * @raymondcast18 made their first contribution in #649.
  * @PhyloevoTi made their first contribution in #655.
  * @ixchelgzlzr made their first contribution in #703.
  * @basanta33 made their first contribution in #704.
  * @Levi-Raskin made their first contribution in #711.

# RevBayes 1.2.5 (Dec 19, 2024)

## Backwards-incompatible changes
  * Remove `underPrior` argument to `mcmc.run( )` and `mcmcmc.run( )`.  You can use `model.ignoreAllData()` instead.
  * Remove `dnLogexponential`. You can use `dnExponential(l) |> tnLog()` instead.
  * Remove `tree.makeBifurcating()` in favor of `tree.suppressOutdegreeOneNodes()`.  To resolve a trifurcation
    at the root you can use `reroot(..., makeBifurcating=TRUE)` instead.

## Features
  * Allow fossil age sampling with an initial tree in `dnConstrainedTopology` (#481).
  * Transformed distributions:
      - Add shift/scale/exp/log/logit/invlogit-transformed distributions (e.g. `dnExp(1) |> tnShift(2)`). (#466, #500)
      - Allow `C + dist`, `C - dist`, `dist - C`, and `C*dist` for a constant `C` and distribution `dist`. (#617)
  * Allow ignoring data from specific nodes with `model.ignoreData(seqData,morphData)` (#628)
  * Extend `type( )` to optionally show the full type-spec (#504).

## Bug fixes
  * Multiple fixes for building, especially on Windows (#490, #491, #493, #494, #495, #496, #498).
  * Checkpointing
      - Fix FBD checkpointing bugs related to `dnBDSTP` (#484) and fossil time moves (#502, #517).
      - Fix MC^3 checkpointing bug (#553).
  * Types
      - Fix type conversion to integer so that it employs deterministic nodes (#545).
      - Make `RealPos` coherent between conversion and construction (#554).
      - Make `vectorFlatten` work on more type (#514).
  * Crash / NaN
      - Fix `dnMixture` of rate matrices and other `Cloneable` objects (#501).
      - Fix `treeAssembly` sometimes failing to initialize branch lengths (#509).
      - Fix crash if a file changes while we are `source( )`-ing it (#510).
      - Fix segfault with `dnConstrainedTopology` + `dnUniformTimeTree` (#513).
      - Fix a segfault in `Model::getOrderedStochasticNodes` (#569).
  * MCMC: initialization
      - Fix fossil age sampling with an initial tree in `dnBDSTP` (#480).
      - Fix FBD initialization issues (#537, #561).
      - Extend the starting tree simulator to account for origin/root age (#561).
  * MCMC
      - Fix recalculation of likelihoods for Brownian motion distributions (#596).
      - Fix recalculation of likelihoods with `mvRootTimeSlideUniform` (#594).
      - Fix recalculation of likelihood with `treeAssembly` (#549).
      - Fix Gelman-Rubin (PSRF) and stationarity stopping rules (#555).
      - Don't move fossil tips outside their age ranges (#559).
      - Fix operator summary for MC^3 when moves are tuned (#522).
  * Debug:
      - Log reason for -Inf and NaN probabilities (#592).
      - Add options `debugMCMC` and `logMCMC` for investigating likelihood recalculation bugs (#570, #575, #593).
      - More informative error message in `treeAssembly` when number of branch lengths is wrong (#605).
  * Misc
      - Allow reading non-square matrices (#564).
      - Allow checking `args.size()` when no arguments are given. (#479).
      - Balance braces when printing `Matrix<Real>` (#615).

## Documentation improvements
  * `dnPhyloCTMC( )` (#487).
  * `fnDiscretizeBeta( )` (#625).
  * `model( )` (#538).
  * `time( )` (#572, #574)
  * `--setOption` / `setOption()` / `getOption()` (#583)
  * Clarify differences between `.clamp()` and `.setValue()` (#599).
  * Stopping and convergence rules (#488).
  * `mcmc` and `mcmcmc`
      - `.run( )` (#485, #488).
      - `.initializeFromCheckpoint( )` (#505).
      - `moveschedule` parameter and the `weight` parameter of moves (#506).
  * Corrections to `dnBivariatePoisson` (#539) and `mcmcmc` (#541).

# RevBayes 1.2.4 (May 29, 2024)

## Features
  * Refactor coalescent to allow heterochronous samples and more.
  * Allow constructing substitution models as SiteMixtureModel objects.
  * Imputation of missing entries in AverageDistanceMatrix

## Bug fixes
  * Fix simulation of rate categories with variable coding.
  * Fix operator[] in AverageDistanceMatrix

# RevBayes 1.2.3 (Apr 26, 2024)

## Changes
  * Setting collapseSampledAncestors=true is now ignored -- use fnCollapseSA(tree).
  * Temporarily disable FBD-Range models (#449).
  * Refactor PhyloOrnsteinUhlenbeckREML (#426)

## Features
  * Allow mnFile and mnModel to write JSON output if given format="json" (#377)
  * Refactor dnBDSTP, dnFBDP, and dnPhylodynamicBDP (#386, #440)
  * Compute average distance matrices more efficiently (#454).
  * Compute ExponentialError probability densities more efficiently (#432, #434).
  * Add member procedure to set node age (#380).
  * Add member procedures to remove invariant and gap sites (#379, #392).
  * Be less picky about initial tree ages matching ages from taxon files (#455).
  * Make mvCollapseExpandFossilBranch complain if r=1 (no fossil sampling) (#440).
  * Make power posterior analyses more flexible (#397).

## Bug fixes
  * Stop BD simulations from hanging in situations with many constraints (#453).
  * Stop moving fossil tip nodes in mvRootTimeSlideUniform (#447).
  * Fix consensusTree (#441).
  * Fix crash reading trees with bad indices (#403, #395).
  * Fix ancestral state reconstruction with missing data (#396)
  * Fix the uniform topology prior (#442).
  * Fix vectorFlatten (#445, #389).
  * Fix a memory leak in reversible jump mixtures (#430).
  * Fix BDSTP segfault (#367)
  * Fix updating phylogenetic likelihood when site mixture probability changes (#437).
  * Fix interaction of multiple stopping rules (#458).
  * Fix crash with tuning on mvScale (#452).
  * Fix writing delimited character data (#362)
  * Fix representation of sampled ancestors in tensorphylo functions (#369).
  * Fix offset for the read tree trace function (#381).
  * Fix Probability(-1) and Probability(2) (#410).
  * Don't flatten arrays of stochastic variables in checkpoint files (#448).


# RevBayes 1.2.2 (Jun 7, 2023)

## Features
  * We now use a flag (-j) instead of a separate executable for Jupyter mode.
  * Matrices can be added, subtracted, and negated.
  * Allow computing SSE likelihoods using TensorPhylo
  * Add pipe operator: `x |> f(y)` now yields `f(x,y)`.
  * ... and many bug fixes.

# RevBayes 1.2.1 "Peitenimi" (Nov 7, 2022)

## Speed & memory
  * Cache transition probability matrices.
  * 4x faster tree summaries.
  * Discrete character data sets take 5x less RAM, and are 5x faster to load.

## Features
  * Add Occurence Birth-Death Process.
  * Better mixing statistics for MC^3.

  * Automatically remove degree-2 root nodes when reading non-clock trees.

## Bug Fixes
  * Fix some problems with the BDSTP.
  * Fix false claim of conflicting constraints. (#250, #288)
  * Don't get stuck on MCMC with amino-acid models.
  * Fix crash in unrooted NNI and SPR.
  * Fix using file paths and directories on Windows.
  * Restore dnFBDP synonym for dnBDSTP.
  * Restore initialTree argument to dnFBDP .
  * Fix check for number of rates in FBD-Range process.
  * Fix crash when reading some trees with sampled ancestors. (#240)
  * Prevent NumStates methods from overriding character exclusion. (#252)
  * C++ standard switched to C++17.

# RevBayes 1.2.0 "Lagerstätte" (Aug 3, 2022)

## New MCMC moves

 * slice sampling (mvSlice) can now do geometric scaling.

## New models/analyses

 * New codon models:
   * Goldman-Yang (1994), Muse-Gaut (1994), FMutSel (Yang 2008)
   * Create of Codon models using stackable blocks fndNdS, fnX3, fnMutSel, fnMutSelAA
   * i.e. fnMutSel( fndNdS( fnX3( fnGTR(er, pi), omega ), F ) # GTR + X3 + dNdS + MutSel
 * New dinucleotide models: fnX2, fnMutSel
   * i.e. fnMutSel( fnX2( fnHKY(kappa, nuc_pi) ), F )    # HKY + X2 + MutSel
 * New birth death models:
   * Birth death sampling treatment process (Magee et al. 2020)
   * Time-heterogeneous fossilized birth death range process
  
## New features

  * optional weighting in fnAverageDistanceMatrix
  * FBD range model uses individual fossil occurrence data
  * site mixture allocations work with codon models
  * automatic handling of whitespace delimited files

## Bug fixes

  * fixed occasional crashes when using recovering a tree from a checkpoint file due to rounding. Checkpointing now records doubles without rounding
  * fixed a bug when printing individual elements of average distance matrices and calculating their completeness


# RevBayes 1.1.1

## Changes

 * **Warning**: this version includes changes to the Random Number Generator, meaning analysis output will be different from v1.1.0, even when run with the same seed.

## New models/analyses

 * new PoMo rate matrices with selection

## Bug fixes

 * fixed an error in the likelihood calculation for discrete characters when using ambiguous characters and partitioning standard data by state number
 * fixed a bug causing unbounded likelihoods when using PhyloBrownianProcessREML with unrooted trees
 * fixed a crash when using ordered rate matrices
 * fixed a crash when printing the value of a rate matrix
 * fixed a crash when using nested clade constraints
 * fixed a crash when using DEC models
 * fixed an issue with using simplex indexes in deterministic nodes
 * fixed a bug allowing NodeTimeSlideUniform to move sampled ancestor nodes
