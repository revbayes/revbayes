/**
 * @file
 * This file contains the Workspace function that adds types and functions
 * to the global workspace, registering them with the interpreter/compiler
 * during the process.
 *
 * @brief Function registering language objects
 *
 * Instructions
 *
 * This is the central registry of Rev objects. It is a large file and needs
 * to be properly organized to facilitate maintenance. Follow these simple
 * guidelines to ensure that your additions follow the existing structure.
 *
 * 1. All headers are added in groups corresponding to directories in the
 *    revlanguage code base.
 * 2. All objects (types, distributions, and functions) are registered in
 *    groups corresponding to directories in the revlanguage code base.
 * 3. All entries in each group are listed in alphabetical order.
 *
 * Some explanation of the directory structure is provided in the comments
 * in this file. Consult these comments if you are uncertain about where
 * to add your objects in the code.
 */

#include <sstream>
#include <vector>
#include <cstdlib>
#include <cstdio>

/* Files including helper classes */
#include "RbException.h"
#include "RlUserInterface.h"
#include "Workspace.h"

/// Miscellaneous types ///

#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
//#include "NumUniqueInVector.h" //suggested by IWYU but breaks the build
#include "RbVector.h"
#include "RevPtr.h"
#include "RlConstantNode.h"
#include "RlDeterministicNode.h"
#include "RlTypedDistribution.h"
#include "RlTypedFunction.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

/* Base types (in folder "datatypes") */

/* Primitive types (in folder "datatypes/basic") */
#include "Integer.h"
#include "IntegerPos.h"
#include "Natural.h"
#include "Probability.h"
#include "Real.h"
#include "RealPos.h"

/* Container types (in folder "datatypes/container") */
#include "ModelVector.h"

/* Taxon types (in folder "datatypes/evolution") */

/* Math types (in folder "datatypes/math") */
#include "RlSimplex.h"

/* Argument rules (in folder "functions/argumentrules") */

/* Basic functions (in folder "functions/basic"). */

/* These are core functions for the Rev environment, providing user help
   and other essential services. */


/* Functions related to evolution (in folder "functions/phylogenetics") */
#include "Func_AlleleFrequencySimulator.h"
#include "Func_AlleleFrequencyMatrixSimulator.h"
#include "Func_avgDistanceMatrix.h"
#include "Func_BirthDeathSimulator.h"
#include "Func_branchScoreDistance.h"
#include "Func_checkNodeOrderConstraints.h"
#include "Func_chronoToPhylo.h"
#include "Func_computeWeightedNodeOrderConstraintsScore.h"
#include "Func_combineCharacter.h"
#include "Func_concatenate.h"
#include "Func_concatenateContinuousCharacterData.h"
#include "Func_CladeSpecificHierarchicalBranchRate.h"
#include "Func_concatenateFromVector.h"
#include "Func_constructRootedTripletDistribution.h"
#include "Func_earlyBurstRates.h"
#include "Func_extantTree.h"
#include "Func_formatDiscreteCharacterData.h"
#include "Func_inferAncestralPopSize.h"
#include "Func_maximumTree.h"
#include "Func_mrcaIndex.h"
#include "Func_nodeAgeByID.h"
#include "Func_phyloDiversity.h"
#include "Func_PhylogeneticIndependentContrasts.h"
#include "Func_PhylogeneticIndependentContrastsMultiSample.h"
#include "Func_pomoStateConverter.h"
#include "Func_pomoRootFrequencies.h"
#include "Func_pruneTree.h"
#include "Func_collapseSA.h"
#include "Func_featureInformedRates.h"
#include "Func_simStartingTree.h"
#include "Func_simTree.h"
#include "Func_simCompleteTree.h"
#include "Func_stitchTree.h"
#include "Func_symmetricDifference.h"
#include "Func_tmrca.h"
#include "Func_treeAssembly.h"
#include "Func_treePairwiseDistances.h"
#include "Func_treePairwiseNodalDistances.h"
#include "Func_treeScale.h"
#include "Func_UPGMA.h"
#include "Func_PoMoStationaryFrequencies.h"
#include "Func_PoMoReversibleMutationRates.h"
#include "Func_FlowT2Populations.h"

/* Frequency functions (in folder "functions/phylogenetics/frequencies") */
#include "Func_F1x4.h"
#include "Func_F3x4.h"
#include "Func_F2x4.h"

/* Rate matrix functions (in folder "functions/phylogenetics/ratematrix") */
#include "Func_BinaryMutationCoalescentRateMatrix.h"
#include "Func_blosum62.h"
#include "Func_biogeographyRateMatrix.h"
#include "Func_chromosomes.h"
#include "Func_chromosomesPloidy.h"

#include "Func_GammaRateModel.h"
#include "Func_InvModel.h"
#include "Func_MixtureModel.h"
#include "Func_ScaleSiteMixtureModel.h"
#include "Func_ScaleVectorSiteMixtureModel.h"
#include "Func_UnitMixture.h"
#include "Func_ConvertRateMatrix.h"
#include "Func_ConvertVectorRateMatrix.h"

#include "Func_codonSynonymousNonsynonymousRateMatrix.h"
#include "Func_codonSynonymousNonsynonymousHKYRateMatrix.h"
#include "Func_GoldmanYang94RateMatrix.h"
#include "Func_MuseGaut94RateMatrix.h"
#include "Func_MuseGaut94KRateMatrix.h"
#include "Func_X3RateMatrix.h"
#include "Func_dNdSRateMatrix.h"
#include "Func_FMutSelRateMatrix.h"
#include "Func_FMutSel0RateMatrix.h"
#include "Func_MutSelRateMatrix.h"
#include "Func_MutSelAARateMatrix.h"

#include "Func_X2RateMatrix.h"

#include "Func_covarionRateMatrix.h"
#include "Func_covarion.h"
#include "Func_cpRev.h"
#include "Func_dayhoff.h"
#include "Func_DECRateMatrix.h"
#include "Func_epoch.h"
#include "Func_f81.h"
#include "Func_FreeBinary.h"
#include "Func_FreeK.h"
#include "Func_freeSymmetricRateMatrix.h"
#include "Func_gtr.h"
#include "Func_hky.h"
#include "Func_hiddenStateRateMatrix.h"
#include "Func_InfiniteSitesRateMatrix.h"
#include "Func_jc.h"
#include "Func_jones.h"
#include "Func_k80.h"
#include "Func_Kimura81.h"
#include "Func_lg.h"
#include "Func_mtRev.h"
#include "Func_mtMam.h"
#include "Func_orderedRateMatrix.h"
#include "Func_PoMoKN.h"
#include "Func_PoMoKNrecurrentMutations.h"
#include "Func_PoMoBalanceKN.h"
#include "Func_revPoMoKN.h"
#include "Func_revPoMoBalanceKN.h"
#include "Func_revPoMoM2N.h"
#include "Func_rtRev.h"
#include "Func_vt.h"
#include "Func_t92.h"
#include "Func_TamuraNei.h"
#include "Func_TIM.h"
#include "Func_TVM.h"
#include "Func_wag.h"

/* Functions related to evolution (in folder "functions/popgen") */
#include "Func_PattersonsD.h"
#include "Func_SegregatingSites.h"
#include "Func_TajimasD.h"
#include "Func_TajimasPi.h"
#include "Func_WattersonTheta.h"


/* Rate map functions (in folder "functions/evolution/ratemap") */
#include "Func_adjacentRateModifier.h"
#include "Func_biogeo_de.h"
#include "Func_distanceRateModifier.h"
#include "Func_generalRateGeneratorSequence.h"
#include "Func_rangeEvolutionRateModifier.h"
#include "Func_stateCountRateModifier.h"
#include "Func_siteRateModifier.h"
#include "Func_hostSwitchRateModifier.h"


/* Cladogeneic state prob function */
#include "Func_biogeographyCladoEventsBD.h"
#include "Func_DECCladoProbs.h"
#include "Func_DECRates.h"
#include "Func_DECRoot.h"
#include "Func_EpochCladoProbs.h"
#include "Func_chromosomesCladoProbs.h"
#include "Func_chromosomesCladoEventsBD.h"
#include "Func_chromosomesPloidyCladoEventsBD.h"
#include "Func_cladogeneticSpeciationRateMatrix.h"
#include "Func_cladogeneticProbabilityMatrix.h"
#include "Func_MixtureCladoProbs.h"
#include "Func_SampledCladogenesisRootFrequencies.h"


/* Input/output functions (in folder "functions/io") */
#include "Func_readPoMoCountFile.h"
#include "Func_convertCountFileToNaturalNumbers.h"
#include "Func_convertFastaFileToNaturalNumbers.h"
#include "Func_convertVCFtoCountsFile.h"

/* Math functions (in folder "functions/math") */
#include "Func_abs.h"
#include "Func_absVector.h"
#include "Func_ceil.h"
#include "Func_choose.h"
#include "Func_coala.h"
#include "Func_diagonalMatrix.h"
#include "Func_empiricalQuantile.h"
#include "Func_exp.h"
#include "Func_expVector.h"
#include "Func_floor.h"
#include "Func_gamma.h"
#include "Func_lnProbability.h"
#include "Func_geographicalDistance.h"
#include "Func_geometricMean.h"
#include "Func_hyperbolicTangent.h"
#include "Func_hyperbolicSine.h"
#include "Func_ln.h"
#include "Func_log.h"
#include "Func_logit.h"
#include "Func_logistic.h"
#include "Func_matrix.h"
#include "Func_max.h"
#include "Func_mean.h"
#include "Func_meanPositive.h"
#include "Func_meanSimplex.h"
#include "Func_median.h"
#include "Func_min.h"
#include "Func_normalize.h"
#include "Func_posteriorPredictiveProbability.h"
//#include "Func_power.h"
//#include "Func_powerVector.h"
#include "Func_round.h"
#include "Func_shortestDistance.h"
#include "Func_sigmoid.h"
#include "Func_sigmoidVector.h"
#include "Func_SmoothTimeLine.h"
#include "Func_sort.h"
#include "Func_sum.h"
#include "Func_sumPositive.h"
#include "Func_sumInteger.h"
#include "Func_sumNatural.h"
#include "Func_standardDeviation.h"
#include "Func_sqrt.h"
#include "Func_trunc.h"
#include "Func_upperTriangle.h"
#include "Func_variance.h"
#include "Func_vectorFlatten.h"


/* Statistics functions (in folder "functions/statistics") */
/* These are functions related to statistical distributions */
#include "Func_assembleContinuousMRF.h"
#include "Func_betaBrokenStick.h"
#include "Func_discretizeBeta.h"
#include "Func_discretizeBetaQuadrature.h"
#include "Func_discretizeGamma.h"
#include "Func_discretizeGammaFromBetaQuantiles.h"
#include "Func_discretizeGammaQuadrature.h"
#include "Func_discretizeLognormalQuadrature.h"
#include "Func_discretizeDistribution.h"
#include "Func_discretizePositiveDistribution.h"
#include "Func_discretizeProbabilityDistribution.h"
#include "Func_dppConcFromMean.h"
#include "Func_dppMeanFromConc.h"
#include "Func_fnNormalizedQuantile.h"
#include "Func_numUniqueInVector.h"
#include "Func_rateShifts.h"
#include "Func_stirling.h"
#include "Func_varianceCovarianceMatrix.h"
#include "Func_decomposedVarianceCovarianceMatrix.h"
#include "Func_partialToCorrelationMatrix.h"

/* Type conversions */
#include "Proc_StringToInt.h"


/** Initialize global workspace */
void RevLanguage::Workspace::initializeFuncGlobalWorkspace(void)
{
    try
    {
        ///////////////////////////////////////////
        /* Add functions (in "functions" folder) */
        ///////////////////////////////////////////

        addFunction( new Func_FlowT2Populations()      );


        /* Rate matrix generator functions (in folder "functions/evolution/ratematrix") */
        addFunction( new Func_BinaryMutationCoalescentRateMatrix()          );
        addFunction( new Func_blosum62()                                    );
        addFunction( new Func_biogeographyRateMatrix()                      );
        addFunction( new Func_chromosomes()                                 );
        addFunction( new Func_chromosomesPloidy()                           );

        addFunction( new Func_ConvertRateMatrix()                           );
        addFunction( new Func_ConvertVectorRateMatrix()                     );

        addFunction( new Func_GammaRateModel()                              );
        addFunction( new Func_InvModel()                                    );
        addFunction( new Func_MixtureModel()                                );
        addFunction( new Func_UnitMixture()                                 );
        addFunction( new Func_ScaleSiteMixtureModel()                       );
        addFunction( new Func_ScaleVectorSiteMixtureModel()                 );

        addFunction( new Func_codonSynonymousNonsynonymousRateMatrix()      );
        addFunction( new Func_codonSynonymousNonsynonymousHKYRateMatrix()   );
        addFunction( new Func_GoldmanYang94RateMatrix()                     );
        addFunction( new Func_MuseGaut94RateMatrix()                        );
        addFunction( new Func_MuseGaut94KRateMatrix()                       );
        addFunction( new Func_X3RateMatrix()                                );
        addFunction( new Func_dNdSRateMatrix()                              );
        addFunction( new Func_FMutSelRateMatrix()                           );
        addFunction( new Func_FMutSel0RateMatrix()                          );
        addFunction( new Func_MutSelRateMatrix()                            );
        addFunction( new Func_MutSelAARateMatrix()                          );

        addFunction( new Func_X2RateMatrix()                                );

        addFunction( new Func_covarionRateMatrix()                          );
        addFunction( new Func_covarion()                                    );
        addFunction( new Func_cpRev()                                       );
        addFunction( new Func_dayhoff()                                     );
        addFunction( new Func_DECRateMatrix()                               );
        addFunction( new Func_epoch()                                       );
        addFunction( new Func_f81()                                         );
        addFunction( new Func_FreeBinary()                                  );
        addFunction( new Func_FreeK()                                       );
        addFunction( new Func_freeSymmetricRateMatrix()                     );
        addFunction( new Func_gtr()                                         );
        addFunction( new Func_hky()                                         );
        addFunction( new Func_hiddenStateRateMatrix()                       );
        addFunction( new Func_InfiniteSitesRateMatrix()                     );
        addFunction( new Func_jc()                                          );
        addFunction( new Func_jones()                                       );
        addFunction( new Func_k80()                                         );
        addFunction( new Func_Kimura81()                                    );
        addFunction( new Func_lg()                                          );
        addFunction( new Func_mtMam()                                       );
        addFunction( new Func_mtRev()                                       );
        addFunction( new Func_orderedRateMatrix()                           );
        addFunction( new Func_PoMoKN()                                      );
        addFunction( new Func_PoMoKNrecurrentMutations()                    );
        addFunction( new Func_PoMoBalanceKN()                               );
        addFunction( new Func_revPoMoKN()                                   );
        addFunction( new Func_revPoMoBalanceKN()                            );
        addFunction( new Func_revPoMoM2N()                                  );
        addFunction( new Func_rtRev()                                       );
        addFunction( new Func_t92()                                         );
        addFunction( new Func_TamuraNei()                                   );
        addFunction( new Func_TIM()                                         );
        addFunction( new Func_TVM()                                         );
        addFunction( new Func_vt()                                          );
        addFunction( new Func_wag()                                         );

        /* frequency functions (in folder "function/phylogenetics/frequencies" */
        addFunction( new Func_F1x4()                                        );
        addFunction( new Func_F3x4()                                        );
        addFunction( new Func_F2x4()                                        );

        /* rate maps used for data augmentation (in folder "functions/evolution/ratemap") */
        addFunction( new Func_adjacentRateModifier() );
        addFunction( new Func_biogeo_de() );
        addFunction( new Func_distanceRateModifier() );
        addFunction( new Func_generalRateGeneratorSequence() );
        addFunction( new Func_hostSwitchRateModifier() );
        addFunction( new Func_rangeEvolutionRateModifier() );
        addFunction( new Func_stateCountRateModifier() );
        addFunction( new Func_siteRateModifier() );

        /* cladogenic probs used for e.g. DEC models (in folder "functions/phylogenetics") */
        addFunction( new Func_avgDistanceMatrix() );
        addFunction( new Func_DECCladoProbs() );
        addFunction( new Func_DECRates() );
        addFunction( new Func_DECRoot() );
        addFunction( new Func_EpochCladoProbs() );
        addFunction( new Func_biogeographyCladoEventsBD() );
        addFunction( new Func_chromosomesCladoProbs() );
        addFunction( new Func_chromosomesCladoEventsBD() );
        addFunction( new Func_chromosomesPloidyCladoEventsBD() );
        addFunction( new Func_CladeSpecificHierarchicalBranchRate() );
        addFunction( new Func_cladogeneticSpeciationRateMatrix() );
        addFunction( new Func_cladogeneticProbabilityMatrix() );
        addFunction( new Func_MixtureCladoProbs() );
        addFunction( new Func_SampledCladogenesisRootFrequencies() );
        addFunction( new Func_PoMoStationaryFrequencies() );
        addFunction( new Func_PoMoReversibleMutationRates() );


		/* Functions related to phylogenetic trees (in folder "functions/phylogenetics/tree") */
        addFunction( new Func_AlleleFrequencySimulator()                        );
        addFunction( new Func_AlleleFrequencyMatrixSimulator()                  );
        addFunction( new Func_BirthDeathSimulator()                             );
        addFunction( new Func_branchScoreDistance()                             );
        addFunction( new Func_checkNodeOrderConstraints()                       );
        addFunction( new Func_chronoToPhylo()                                   );
        addFunction( new Func_computeWeightedNodeOrderConstraintsScore()        );
        addFunction( new Func_combineCharacter()                                );
        addFunction( new Func_concatenate()                                     );
        addFunction( new Func_concatenateContinuousCharacterData()              );
        addFunction( new Func_concatenateFromVector()                           );
        addFunction( new Func_constructRootedTripletDistribution()              );
        addFunction( new Func_formatDiscreteCharacterData()                     );
        addFunction( new Func_EarlyBurstRates()                                 );
        addFunction( new Func_extantTree()                                      );
        addFunction( new Func_inferAncestralPopSize()                           );
        addFunction( new Func_maximumTree()                                     );
        addFunction( new Func_mrcaIndex()                                       );
        addFunction( new Func_nodeAgeByID()                                     );
        addFunction( new Func_phyloDiversity()                                  );
        addFunction( new Func_PhylogeneticIndependentContrasts()                );
        addFunction( new Func_PhylogeneticIndependentContrastsMultiSample()     );
        addFunction( new Func_pomoStateConverter()                              );
        addFunction( new Func_pomoRootFrequencies()                             );
        addFunction( new Func_pruneTree()                                       );
        addFunction( new Func_collapseSA<Tree>()                                );
        addFunction( new Func_collapseSA<BranchLengthTree>()                    );
        addFunction( new Func_collapseSA<TimeTree>()                            );
        addFunction( new Func_featureInformedRates()                            );
        addFunction( new Func_readPoMoCountFile()                               );
        addFunction( new Func_convertCountFileToNaturalNumbers()                );
        addFunction( new Func_convertFastaFileToNaturalNumbers()                );
        addFunction( new Func_convertVCFtoCountsFile()                          );
        addFunction( new Func_simStartingTree()                                 );
        addFunction( new Func_simTree()                                         );
        addFunction( new Func_simCompleteTree()                                 );
        addFunction( new Func_stitchTree()                                      );
        addFunction( new Func_symmetricDifference()                             );
        addFunction( new Func_tmrca()                                           );
        addFunction( new Func_treePairwiseDistances()                           );
        addFunction( new Func_treePairwiseNodalDistances()                      );
        addFunction( new Func_treeAssembly()                                    );
        addFunction( new Func_treeScale()                                       );
        addFunction( new Func_UPGMA()                                           );

        /* Population genetics functions (in folder "functions/popgen") */
        addFunction( new Func_PattersonsD()      );
        addFunction( new Func_SegregatingSites() );
        addFunction( new Func_TajimasD()         );
        addFunction( new Func_TajimasPi()        );
        addFunction( new Func_WattersonTheta()   );


        /* Math functions (in folder "functions/math") */

		// absolute function
        addFunction( new Func_abs()                  );
        addFunction( new Func_absVector()            );

		// ceil function
        addFunction( new Func_ceil<Real,Integer>()  );
        addFunction( new Func_ceil<RealPos,Natural>()  );

        // choose function
        addFunction( new Func_choose() );

        // coala function
        addFunction( new Func_coala()        );

        // diagonal matrix
        addFunction( new Func_diagonalMatrix() );

        // empirical quantile function
        addFunction( new Func_empiricalQuantile()  );

        // exponential function
        addFunction( new Func_exp() );
        addFunction( new Func_expVector() );

		// floor function
        addFunction( new Func_floor<Real,Integer>()  );
        addFunction( new Func_floor<RealPos,Natural>()  );

        // gamma function
        addFunction( new Func_gamma() );

        // geometric mean function
        addFunction( new Func_geometricMean() );

        // logistic function
        addFunction( new Func_logistic() );

        // natural log function
        addFunction( new Func_ln()  );

        // log function
        addFunction( new Func_log()  );

        // logit function
        addFunction( new Func_logit()  );

        // matrix function (converts into MatrixReal)
        addFunction( new Func_matrix() );

        // max function
        addFunction( new Func_max()  );

        // mean function
        addFunction( new Func_mean()  );
        addFunction( new Func_meanPositive()  );
        addFunction( new Func_meanSimplex()  );

        // median function
        addFunction( new Func_median()  );

        // min function
		addFunction( new Func_min()  );

        // normalize vector function
		addFunction( new Func_normalize()  );

        // round function
        addFunction( new Func_round<Real,Integer>()  );
        addFunction( new Func_round<RealPos,Natural>()  );

        // sort vector function
        addFunction( new Func_sort() );

        // sigmoid function
        addFunction( new Func_sigmoid() );
        addFunction( new Func_sigmoidVector() );

        // rate shift function
        addFunction( new Func_shiftEvents<RealPos>()              );
        addFunction( new Func_shiftEvents<ModelVector<RealPos>>() );

		// square root function
        addFunction( new Func_sqrt()  );

        // sum function
        addFunction( new Func_sum()  );
        addFunction( new Func_sumPositive()  );
        addFunction( new Func_sumInteger()  );
        addFunction( new Func_sumNatural()  );

        // standard deviation function
        addFunction( new Func_standardDeviation()  );

        // geographical distance function
        addFunction( new Func_geographicalDistance() );
        addFunction( new Func_shortestDistance() );

                // hyperbolic tangent function
        addFunction( new Func_hyperbolicTangent() );

        // hyperbolic sine function
        addFunction( new Func_hyperbolicSine() );

		// truncate function
        addFunction( new Func_trunc<Real,Integer>()  );
        addFunction( new Func_trunc<RealPos,Natural>()  );

        // upper triangle of a matrix function
        addFunction( new Func_upperTriangle()  );

        // variance function
        addFunction( new Func_variance()  );

        // vector flatten
        addFunction( new Func_vectorFlatten<Real>() );
        addFunction( new Func_vectorFlatten<RealPos>() );
        addFunction( new Func_vectorFlatten<Probability>() );
        addFunction( new Func_vectorFlatten<Integer>() );
        addFunction( new Func_vectorFlatten<Natural>() );

        // get ln Probability function
        addFunction( new Func_lnProbability() );

        // empirical cummulative probability function
        addFunction( new Func_posteriorPredictiveProbability()  );


        /* Statistics functions (in folder "functions/statistics") */

        // helpers for Markov Random Field models
        addFunction( new Func_assembleContinuousMRF( )     );
        addFunction( new Func_SmoothTimeLine( )     );

		// some helper statistics for the DPP distribution
        addFunction( new Func_dppConcFromMean( )     );
        addFunction( new Func_dppMeanFromConc( )  );
        addFunction( new Func_stirling( )     );

		// count the number of unique elements in vector
        addFunction( new Func_numUniqueInVector<Real>( )  );
        addFunction( new Func_numUniqueInVector<RealPos>( )  );
        addFunction( new Func_numUniqueInVector<Integer>( )  );
        addFunction( new Func_numUniqueInVector<Natural>( )  );
        addFunction( new Func_numUniqueInVector<Probability>( )  );
        addFunction( new Func_numUniqueInVector<Simplex>( )  );

        // return a distcretized (by quantile) and normalized vector from a continuous distribution
        addFunction( new Func_fnNormalizedQuantile<Real>()    );
        addFunction( new Func_fnNormalizedQuantile<RealPos>()    );
        
        addFunction( new Func_discretizeDistribution( )            );
        addFunction( new Func_discretizePositiveDistribution( )    );
        addFunction( new Func_discretizeProbabilityDistribution( ) );

        // return a discretized gamma distribution (for gamma-dist rates)
        addFunction( new Func_discretizeBeta( )    );
        addFunction( new Func_discretizeBetaQuadrature( )    );
        addFunction( new Func_discretizeGamma( )   );
        addFunction( new Func_discretizeGammaFromBetaQuantiles( )   );
        addFunction( new Func_discretizeGammaQuadrature( )   );
        addFunction( new Func_discretizeLognormalQuadrature( )   );

        addFunction( new Func_betaBrokenStick( )   );

        addFunction( new Func_varianceCovarianceMatrix( )           );
        addFunction( new Func_decomposedVarianceCovarianceMatrix( ) );
        addFunction( new Func_partialToCorrelationMatrix( )         );


        // Type conversion
        addFunction( new Proc_StringToInt( )                         );

    }
    catch(RbException& rbException)
    {

        RBOUT("Caught an exception while initializing functions in the workspace\n");
        std::ostringstream msg;
        rbException.print(msg);
        msg << std::endl;
        RBOUT(msg.str());

        RBOUT("Please report this bug to the RevBayes Development Core Team");

        RBOUT("Press any character to exit the program.");
        getchar();
        exit(1);
    }

}
