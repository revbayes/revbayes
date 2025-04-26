# RevBayes 1.3.0 (unreleased)
## What's Changed
* simmap bugfix for when some characters are excluded by @kopperud in https://github.com/revbayes/revbayes/pull/636
* Fix `.dropTip()` behavior for multifurcations by @davidcerny in https://github.com/revbayes/revbayes/pull/640
* Print convergence statistics during MCMC runs with stopping rules by @davidcerny in https://github.com/revbayes/revbayes/pull/645
* Random resolution of multifurcations and MBL time-scaling by @davidcerny in https://github.com/revbayes/revbayes/pull/641
* Add cosine function by @bredelings in https://github.com/revbayes/revbayes/pull/647
* MPI checkpointing fix by @davidcerny in https://github.com/revbayes/revbayes/pull/646
* Prevent loss of integer precision in MC^3 by @davidcerny in https://github.com/revbayes/revbayes/pull/652
* Improvements to stopping rule messages by @davidcerny in https://github.com/revbayes/revbayes/pull/657
* Use different seeds for replicates in MPI run by @bjoelle in https://github.com/revbayes/revbayes/pull/662
* Prevent modification of clamped value by @ms609 in https://github.com/revbayes/revbayes/pull/600
* Prevent NaN heats in MC^3 by @davidcerny in https://github.com/revbayes/revbayes/pull/661
* Dev hoehna lab by @hoehna in https://github.com/revbayes/revbayes/pull/682
* Additional coalescent development by @hoehna in https://github.com/revbayes/revbayes/pull/473
* Don't treat deterministic nodes as constants by @bredelings in https://github.com/revbayes/revbayes/pull/678
* Use std::shared_ptr<Environment> to hold Environments. by @bredelings in https://github.com/revbayes/revbayes/pull/687
* Fix state space when subsetting by state space by @hoehna in https://github.com/revbayes/revbayes/pull/696
* Better documentation and progress monitoring for power posterior analyses by @davidcerny in https://github.com/revbayes/revbayes/pull/673
* Added sine function by @sigibrock in https://github.com/revbayes/revbayes/pull/648
* Fix mvRateAgeBetaShift by @bredelings in https://github.com/revbayes/revbayes/pull/688
* Fix windows long by @bredelings in https://github.com/revbayes/revbayes/pull/708
* BUG #668 dnUniformTopology ignores neg clade const by @Levi-Raskin in https://github.com/revbayes/revbayes/pull/711
* Fix out of range proposal by @bredelings in https://github.com/revbayes/revbayes/pull/709
* Fix deterministic index by @bredelings in https://github.com/revbayes/revbayes/pull/719
* Only type-convert variables based on value if they are constant by @bredelings in https://github.com/revbayes/revbayes/pull/721
* Fix logging when running Mcmcmc twice. by @bredelings in https://github.com/revbayes/revbayes/pull/727
* MC^3 checkpointing fix by @davidcerny in https://github.com/revbayes/revbayes/pull/670
* Fix `readTrees(text = ...)` by @davidcerny in https://github.com/revbayes/revbayes/pull/735
* Add an `nruns` argument to `readTrace()` by @davidcerny in https://github.com/revbayes/revbayes/pull/736
* Fix caching bugs in phyloCTMC by @bredelings in https://github.com/revbayes/revbayes/pull/729
* Allow MCMC moves with low prob ratio and high Jacobian/Hastings ratio by @bredelings in https://github.com/revbayes/revbayes/pull/728
* Display current stopping rule values when resuming from a checkpoint by @davidcerny in https://github.com/revbayes/revbayes/pull/739

## Documentation
* Filled out help file for exp function. by @brpetrucci in https://github.com/revbayes/revbayes/pull/653
* Document simplex moves by @ms609 in https://github.com/revbayes/revbayes/pull/606
* Changes to rate matrix model help files by @brpetrucci in https://github.com/revbayes/revbayes/pull/656
* Completing and standardizing stats distributions help files, part 1 by @brpetrucci in https://github.com/revbayes/revbayes/pull/663
* Document DPP moves by @davidcerny in https://github.com/revbayes/revbayes/pull/666
* add help for jukes cantor by @raymondcast18 in https://github.com/revbayes/revbayes/pull/649
* Created matrix.md and var.md files for review by @raymondcast18 in https://github.com/revbayes/revbayes/pull/681
* Added function covarion documentation by @basanta33 in https://github.com/revbayes/revbayes/pull/704
* Adding documentation to 'fnDiscretizeDistribution' and 'fnDiscretizeG… by @basanta33 in https://github.com/revbayes/revbayes/pull/707
* updated fnFreeBinary help file by @sigibrock in https://github.com/revbayes/revbayes/pull/706
* Wrote help files for power and write function by @raymondcast18 in https://github.com/revbayes/revbayes/pull/690
* vectorFlatten.md fnFreeK.md, fnReadVCF.md, dnUniformInteger.md help files by @raymondcast18 in https://github.com/revbayes/revbayes/pull/667
* Updating documentation for readTreeTrace by @ixchelgzlzr in https://github.com/revbayes/revbayes/pull/703
* Editing WAG help file  by @PhyloevoTi in https://github.com/revbayes/revbayes/pull/655
* Update readTrace.md by @bredelings in https://github.com/revbayes/revbayes/pull/702
* Document TreeTrace methods by @ms609 in https://github.com/revbayes/revbayes/pull/722
* Document readTrace methods by @ms609 in https://github.com/revbayes/revbayes/pull/723

## Infrastructure
* Fix CI builds by dropping openlibm by @bredelings in https://github.com/revbayes/revbayes/pull/644
* Changes to build files to streamline help infrastructure by @brpetrucci in https://github.com/revbayes/revbayes/pull/659
* Revamped tutorial testing by @brpetrucci in https://github.com/revbayes/revbayes/pull/674
* Ensuring integration tests can run tutorial checkpoint tests by @brpetrucci in https://github.com/revbayes/revbayes/pull/676
* Changed website submodule to pull from source rather than master. by @brpetrucci in https://github.com/revbayes/revbayes/pull/692

## Ignore
* Bump version after v1.2.5 release by @bredelings in https://github.com/revbayes/revbayes/pull/637
* Revert accidentally committed-by-unreview patch by @bredelings in https://github.com/revbayes/revbayes/pull/654
* Remove empty touch( ), keep( ), and  restore( ) specializations. by @bredelings in https://github.com/revbayes/revbayes/pull/658
* Update nlohmann::json to avoid tons of compilation warnings with GCC 15. by @bredelings in https://github.com/revbayes/revbayes/pull/685
* Fix broken overrides part1 by @bredelings in https://github.com/revbayes/revbayes/pull/684
* Mark overridden virtual functions as such in tree event moves by @davidcerny in https://github.com/revbayes/revbayes/pull/669
* Don't throw exception to perform a simple check. by @bredelings in https://github.com/revbayes/revbayes/pull/677

## New Contributors
* @raymondcast18 made their first contribution in https://github.com/revbayes/revbayes/pull/649
* @sigibrock made their first contribution in https://github.com/revbayes/revbayes/pull/648
* @basanta33 made their first contribution in https://github.com/revbayes/revbayes/pull/704
* @ixchelgzlzr made their first contribution in https://github.com/revbayes/revbayes/pull/703
* @Levi-Raskin made their first contribution in https://github.com/revbayes/revbayes/pull/711
* @PhyloevoTi made their first contribution in https://github.com/revbayes/revbayes/pull/655

# RevBayes 1.2.5 (Dec 19, 2025)

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
      - Balance braces when printing Matrix<Real> (#615).

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
