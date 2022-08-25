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
#include <cstdio>
#include <cstdlib>

/* Files including helper classes */
#include "RbException.h"
#include "RlUserInterface.h"
#include "Workspace.h"

#include "DPPAllocateAuxGibbsMove.h"
#include "DPPTableValueUpdate.h"
#include "GibbsMixtureAllocationProposal.h"
#include "Integer.h"
#include "IntegerPos.h"
#include "MixtureAllocationProposal.h"
#include "ModelObject.h"
#include "ModelVector.h"
#include "Natural.h"
#include "RbVector.h"
#include "Real.h"
#include "RealPos.h"
#include "RevPtr.h"
#include "ReversibleJumpMixtureProposal.h"
#include "RlTree.h"
#include "SimpleProposal.h"
#include "Simplex.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "UPPAllocationProposal.h"
#include "WorkspaceToCoreWrapperObject.h"


/// Types ///

/* These types are needed as template types for the moves */
#include "RlRateGenerator.h"
#include "Probability.h"
#include "RlSimplex.h"


/// Moves ///

/* Move types (in folder "datatypes/inference/moves") (grouped by parameter type) */

/* Moves on real values */
#include "Move_Mirror.h"
#include "Move_MirrorMultiplier.h"
#include "Move_HSRFHyperpriorsGibbs.h"
#include "Move_HSRFIntervalSwap.h"
#include "Move_HSRFUnevenGridHyperpriorsGibbs.h"
#include "Move_SliceSampling.h"
#include "Move_Scale.h"
#include "Move_GammaScale.h"
#include "Move_RandomDive.h"
#include "Move_ScaleBactrian.h"
#include "Move_ScaleBactrianCauchy.h"
#include "Move_Slide.h"
#include "Move_SlideBactrian.h"

/* Moves on probability values */
#include "Move_BetaProbability.h"

/* Compound Moves on Real Values */
#include "Move_AVMVN.h"
#include "Move_UpDownSlide.h"
#include "Move_UpDownSlideBactrian.h"
#include "Move_UpDownTreeScale.h"
#include "Move_LevyJumpSum.h"
#include "Move_LevyJump.h"

/* Moves on integer values */
#include "Move_RandomGeometricWalk.h"
#include "Move_RandomIntegerWalk.h"
#include "Move_RandomNaturalWalk.h"
#include "Move_BinarySwitch.h"


/* Moves on simplices */
#include "Move_DirichletSimplex.h"
#include "Move_BetaSimplex.h"
#include "Move_ElementSwapSimplex.h"

/* Moves on vectors */
#include "Move_ElementScale.h"
#include "Move_ElementSlide.h"
#include "Move_GMRFHyperpriorGibbs.h"
#include "Move_GMRFUnevenGridHyperpriorGibbs.h"
#include "Move_MultipleElementScale.h"
#include "Move_ShrinkExpand.h"
#include "Move_ShrinkExpandScale.h"
#include "Move_SingleElementScale.h"
#include "Move_SingleElementSlide.h"
#include "Move_EllipticalSliceSamplingSimple.h"
#include "Move_SynchronizedVectorFixedSingleElementSlide.h"
#include "Move_VectorBinarySwitch.h"
#include "Move_VectorSingleElementScale.h"
#include "Move_VectorSingleElementSlide.h"
#include "Move_VectorElementSwap.h"
#include "Move_VectorFixedSingleElementSlide.h"
#include "Move_VectorScale.h"
#include "Move_VectorSlide.h"
#include "Move_VectorSlideRecenter.h"

/* Moves on real valued matrices */
#include "Move_MatrixSingleElementScale.h"
#include "Move_MatrixSingleElementSlide.h"

/* Moves on correlation matrices */
#include "Move_ConjugateInverseWishart.h"
#include "Move_CorrelationMatrixUpdate.h"
#include "Move_CorrelationMatrixElementSwap.h"
#include "Move_CorrelationMatrixRandomWalk.h"
#include "Move_CorrelationMatrixSingleElementBeta.h"
#include "Move_CorrelationMatrixSpecificElementBeta.h"
//#include "Move_CorrelationMatrixPartialSingleElementBeta.h"
//#include "Move_CorrelationMatrixReparameterization.h"
//#include "Move_CorrelationMatrixExpansion.h"

/* Moves on random adjacency graphs */
#include "Move_GraphFlipEdge.h"
#include "Move_GraphFlipClique.h"
#include "Move_GraphShiftEdge.h"

/* Moves on continuous character data (real valued matrices) */
#include "Move_ContinuousCharacterDataSlide.h"

/* Moves on discrete character data matrices */
#include "Move_HomeologPhase.h"

///* Moves on covariance matrices */
#include "Move_MatrixRealSymmetricSlide.h"


/* Moves on mixtures (in folder "datatypes/inference/moves/mixture") */
#include "Move_DPPAllocateAuxGibbsMove.h"
#include "Move_DPPGibbsConcentration.h"
#include "Move_DPPTableValueUpdate.h"
#include "Move_MixtureAllocation.h"
#include "Move_GibbsMixtureAllocation.h"
#include "Move_ReversibleJumpSwitchMove.h"
#include "Move_UPPAllocation.h"

// moves for the DPP table values
#include "ScaleProposal.h"
#include "SlideProposal.h"
#include "BetaSimplexProposal.h"


/* Moves on character histories/data augmentation */
// #include "Move_NodeCharacterHistoryRejectionSample.h"
// #include "Move_PathCharacterHistoryRejectionSample.h"
#include "Move_CharacterHistory.h"

#include "Move_ResampleFBD.h"

#include "Move_BirthDeathEventContinuous.h"
#include "Move_BirthDeathEventDiscrete.h"
#include "Move_BirthDeathFromAgeEvent.h"
#include "Move_ContinuousEventScale.h"
#include "Move_DiscreteEventCategoryRandomWalk.h"
#include "Move_EventTimeBeta.h"
#include "Move_EventTimeSlide.h"
#include "Move_GibbsDrawCharacterHistory.h"

#include "Move_MultiValueEventBirthDeath.h"
#include "Move_MultiValueEventScale.h"
#include "Move_MultiValueEventSlide.h"

#include "Move_OrderedEventTimeSlide.h"
#include "Move_OrderedEventBirthDeath.h"
#include "Move_OrderedEventSlide.h"
#include "Move_OrderedEventScale.h"
#include "Move_OrderedEventVectorSlide.h"
#include "Move_OrderedEventVectorScale.h"


/* Moves on continuous phyloprocesses (Brownian, multivariate Brownian, etc) */

/* Tree proposals (in folder "datatypes/inference/moves/tree") */
#include "Move_AddRemoveTip.h"
#include "Move_BurstEvent.h"
#include "Move_BranchLengthScale.h"
#include "Move_CollapseExpandFossilBranch.h"
#include "Move_IndependentTopology.h"
#include "Move_EmpiricalTree.h"
#include "Move_FNPR.h"
#include "Move_FossilTipTimeUniform.h"
#include "Move_FossilTipTimeSlideUniform.h"
#include "Move_GibbsPruneAndRegraft.h"
#include "Move_GibbsSubtreePruneAndRegraft.h"
#include "Move_GibbsSubtreeSwap.h"
#include "Move_LayeredScaleProposal.h"
#include "Move_NarrowExchange.h"
#include "Move_NNIClock.h"
#include "Move_NNINonclock.h"
#include "Move_NodeRateTimeSlideUniform.h"
#include "Move_NodeTimeScale.h"
#include "Move_NodeTimeSlideUniform.h"
#include "Move_NodeTimeSlideUniformAgeConstrained.h"
#include "Move_NodeTimeSlidePathTruncatedNormal.h"
#include "Move_NodeTimeSlideBeta.h"
#include "Move_RateAgeBetaShift.h"
#include "Move_RateAgeProposal.h"
#include "Move_RateAgeSubtreeProposal.h"
#include "Move_RootTimeScaleBactrian.h"
#include "Move_RootTimeSlide.h"
#include "Move_RootTimeSlideUniform.h"
#include "Move_SpeciesNarrowExchange.h"
#include "Move_SpeciesNodeTimeSlideUniform.h"
#include "Move_SpeciesSubtreeScale.h"
#include "Move_SpeciesSubtreeScaleBeta.h"
#include "Move_SpeciesTreeScale.h"
#include "Move_SubtreeScale.h"
#include "Move_SubtreeSwap.h"
#include "Move_SPRNonclock.h"
#include "Move_TreeScale.h"
//#include "Move_WeightedNodeTimeSlide.h"
#include "Move_WeightedBranchLengthScale.h"


#include "Move_NarrowExchangeRateMatrix.h"

#include "Move_IndependentPriorSampler.h"

/** Initialize global workspace */
void RevLanguage::Workspace::initializeMoveGlobalWorkspace(void)
{

    try
    {
        ///////////////////////////////////////////////////////////////////////////////////
        /* Add moves (in folder "datatypes/inference/moves") (grouped by parameter type) */
        ///////////////////////////////////////////////////////////////////////////////////

        /* Regular moves (in folder "datatypes/inference/moves") (grouped by parameter type) */

        /* Moves on real values */
        addType( new Move_Mirror() );
        addType( new Move_MirrorMultiplier() );
        addType( new Move_Scale() );
        addType( new Move_GammaScale() );
        addType( new Move_RandomDive() );
        addType( new Move_ScaleBactrian() );
        addType( new Move_ScaleBactrianCauchy() );
        addType( new Move_Slide() );
        addType( new Move_SlideBactrian() );
        addType( new Move_SliceSampling() );
        addType( new Move_GMRFHyperpriorGibbs() );
        addType( new Move_GMRFUnevenGridHyperpriorGibbs() );
        /* Moves on probability */
        addType( new Move_BetaProbability() );

        /* compound moves */
//        addType("mvUpDownScale",         new Move_UpDownScale() );
        addType( new Move_AVMVN() );
        addType( new Move_UpDownTreeScale() );
        addType( new Move_UpDownSlide() );
        addType( new Move_UpDownSlideBactrian() );

        // compound moves on real values
        addType( new Move_LevyJumpSum() );
        addType( new Move_LevyJump() );

        /* Moves on integer values */
        addType( new Move_RandomIntegerWalk() );
        addType( new Move_RandomNaturalWalk() );
        addType( new Move_RandomGeometricWalk() );
        addType( new Move_BinarySwitch() );

        /* Moves on simplices */
        addType( new Move_DirichletSimplex() );
        addType( new Move_BetaSimplex() );
        addType( new Move_ElementSwapSimplex() );

        /* Moves on vectors */
        addType( new Move_MultipleElementScale() );
        addType( new Move_HSRFHyperpriorsGibbs() );
        addType( new Move_HSRFIntervalSwap() );
        addType( new Move_HSRFUnevenGridHyperpriorsGibbs() );
        addType( new Move_SingleElementSlide() );
        addType( new Move_SingleElementScale() );
        addType( new Move_ShrinkExpand() );
        addType( new Move_ShrinkExpandScale() );
        addType( new Move_VectorBinarySwitch() );
        addType( new Move_VectorScale() );
        addType( new Move_VectorSlide() );
        addType( new Move_VectorSlideRecenter() );
        addType( new Move_ElementScale() );
        addType( new Move_ElementSlide() );
        addType( new Move_VectorSingleElementScale() );
        addType( new Move_VectorSingleElementSlide() );
        addType( new Move_VectorElementSwap<Real>( ) );
        addType( new Move_VectorElementSwap<RealPos>( ) );
        addType( new Move_VectorElementSwap<Natural>( ) );
        addType( new Move_VectorElementSwap<Integer>( ) );
        addType( new Move_VectorElementSwap<Probability>( ) );
        addType( new Move_VectorElementSwap<Simplex>( ) );
        addType( new Move_VectorElementSwap<ModelVector<RealPos> >( ) );
        addType( new Move_VectorElementSwap<ModelVector<Real> >( ) );
        addType( new Move_VectorElementSwap<RateGenerator>( ) );
        addType( new Move_VectorElementSwap<Tree>( ) );
        addType( new Move_VectorFixedSingleElementSlide() );
        addType( new Move_EllipticalSliceSamplingSimple() );

        addType( new Move_SynchronizedVectorFixedSingleElementSlide() );

        /* Moves on matrices of real values */
        addType( new Move_MatrixSingleElementScale() );
        addType( new Move_MatrixSingleElementSlide() );

        /* Moves on matrices of correlations */
        addType( new Move_CorrelationMatrixUpdate()                   );
        addType( new Move_CorrelationMatrixRandomWalk()               );
        addType( new Move_CorrelationMatrixSingleElementBeta()        );
        addType( new Move_CorrelationMatrixSpecificElementBeta()      );
        addType( new Move_CorrelationMatrixElementSwap()              );

//        addType( new Move_CorrelationMatrixPartialSingleElementBeta() );
//        addType( new Move_CorrelationMatrixReparameterization()       );
//        addType( new Move_CorrelationMatrixExpansion()                );

        /* Moves on matrices of real values */
        addType( new Move_MatrixRealSymmetricSlide() );

        /* Moves on matrices of real values */
        addType( new Move_ConjugateInverseWishart() );

        /* Moves of random adjacency graphs */
        addType( new Move_GraphFlipEdge() );
        addType( new Move_GraphFlipClique() );
        addType( new Move_GraphShiftEdge() );

        /* Moves on continuous character data (matrices of real values) */
        addType( new Move_ContinuousCharacterDataSlide() );

        /* Moves on discrete character data matrices */
        addType( new Move_HomeologPhase() );

        /* Moves on mixtures (in folder "datatypes/inference/moves/mixture") */
        addType( new Move_DPPTableValueUpdate<Real>(    new RevBayesCore::SlideProposal( NULL, 1.0 ) ) );
        addType( new Move_DPPTableValueUpdate<RealPos>( new RevBayesCore::ScaleProposal( NULL, 1.0 ) ) );
        addType( new Move_DPPTableValueUpdate<Simplex>( new RevBayesCore::BetaSimplexProposal( NULL, 10.0 ) ) );

//        addType("mvDPPScaleCatVals",                new Move_DPPScaleCatValsMove() );
//        addType("mvDPPScaleCatAllocateAux",         new Move_DPPScaleCatAllocateAux() );
        addType( new Move_DPPAllocateAuxGibbsMove<Real>() );
        addType( new Move_DPPAllocateAuxGibbsMove<RealPos>() );
        addType( new Move_DPPAllocateAuxGibbsMove<Probability>() );
        addType( new Move_DPPAllocateAuxGibbsMove<Integer>() );
        addType( new Move_DPPAllocateAuxGibbsMove<Natural>() );
        addType( new Move_DPPAllocateAuxGibbsMove<Simplex>() );
        addType( new Move_DPPGibbsConcentration( ) );
        addType( new Move_MixtureAllocation<Real>( ) );
        addType( new Move_MixtureAllocation<RealPos>( ) );
        addType( new Move_MixtureAllocation<Natural>( ) );
        addType( new Move_MixtureAllocation<Integer>( ) );
        addType( new Move_MixtureAllocation<Probability>( ) );
        addType( new Move_MixtureAllocation<Simplex>( ) );
        addType( new Move_MixtureAllocation<ModelVector<RealPos> >( ) );
        addType( new Move_MixtureAllocation<ModelVector<Real> >( ) );
        addType( new Move_MixtureAllocation<RateGenerator>( ) );
        addType( new Move_MixtureAllocation<Tree>( ) );
        addType( new Move_GibbsMixtureAllocation<Real>( ) );
        addType( new Move_GibbsMixtureAllocation<RealPos>( ) );
        addType( new Move_GibbsMixtureAllocation<Natural>( ) );
        addType( new Move_GibbsMixtureAllocation<Integer>( ) );
        addType( new Move_GibbsMixtureAllocation<Probability>( ) );
        addType( new Move_GibbsMixtureAllocation<Simplex>( ) );
        addType( new Move_GibbsMixtureAllocation<RateGenerator>( ) );
        addType( new Move_UPPAllocation<RealPos>() );

        addType( new Move_ReversibleJumpSwitch<Real>( )                  );
        addType( new Move_ReversibleJumpSwitch<RealPos>( )               );
        addType( new Move_ReversibleJumpSwitch<Natural>( )               );
        addType( new Move_ReversibleJumpSwitch<Integer>( )               );
        addType( new Move_ReversibleJumpSwitch<Probability>( )           );
        addType( new Move_ReversibleJumpSwitch<Simplex>( )               );
        addType( new Move_ReversibleJumpSwitch<ModelVector<Natural> >( ) );
        addType( new Move_ReversibleJumpSwitch<Tree>( )                  );

        addType( new Move_MultiValueEventBirthDeath()                    );
        addType( new Move_MultiValueEventScale()                         );
        addType( new Move_MultiValueEventSlide()                         );


        addType( new Move_OrderedEventTimeSlide()                        );
        addType( new Move_OrderedEventBirthDeath()                       );
        addType( new Move_OrderedEventSlide<Real>()                      );
        addType( new Move_OrderedEventSlide<RealPos>()                   );
        addType( new Move_OrderedEventSlide<Probability>()               );
        addType( new Move_OrderedEventScale<RealPos>()                   );
        addType( new Move_OrderedEventVectorSlide<ModelVector<Real>>()        );
        addType( new Move_OrderedEventVectorSlide<ModelVector<RealPos>>()     );
        addType( new Move_OrderedEventVectorSlide<ModelVector<Probability>>() );
        addType( new Move_OrderedEventVectorScale<ModelVector<RealPos>>()     );


        addType( new Move_BirthDeathEventContinuous()                    );
        addType( new Move_BirthDeathEventDiscrete()                      );
        addType( new Move_ContinuousEventScale()                         );
        addType( new Move_DiscreteEventCategoryRandomWalk()              );
        addType( new Move_EventTimeBeta()                                );
        addType( new Move_EventTimeSlide()                               );
        addType( new Move_BirthDeathFromAgeEvent()                       );
        addType( new Move_GibbsDrawCharacterHistory()                    );

        /* Tree proposals (in folder "datatypes/inference/moves/tree") */
        addType( new Move_AddRemoveTip()                     );
        addType( new Move_BurstEvent()                       );
        addType( new Move_BranchLengthScale()                );
        addType( new Move_CollapseExpandFossilBranch()       );
        addType( new Move_IndependentTopology()              );
        addType( new Move_EmpiricalTree()                    );
        addType( new Move_FNPR()                             );
        addType( new Move_FossilTipTimeUniform()             );
        addType( new Move_FossilTipTimeSlideUniform()        );
        addType( new Move_GibbsPruneAndRegraft()             );
        addType( new Move_GibbsSubtreePruneAndRegraft()      );
        addType( new Move_GibbsSubtreeSwap()                 );
        addType( new Move_LayeredScaleProposal()             );
        addType( new Move_NarrowExchange()                   );
        addType( new Move_NarrowExchangeRateMatrix()         );
        addType( new Move_NNIClock()                         );
        addType( new Move_NNINonclock()                      );
        addType( new Move_NodeRateTimeSlideUniform()         );
        addType( new Move_NodeTimeScale()                    );
        addType( new Move_NodeTimeSlidePathTruncatedNormal() );
        addType( new Move_NodeTimeSlideUniform()             );
        addType( new Move_NodeTimeSlideUniformAgeConstrained());
        addType( new Move_NodeTimeSlideBeta()                );
        addType( new Move_RateAgeBetaShift()                 );
        addType( new Move_RateAgeProposal()                  );
        addType( new Move_RateAgeSubtreeProposal()           );
        addType( new Move_RootTimeSlide()                    );
        addType( new Move_RootTimeScaleBactrian()            );
        addType( new Move_RootTimeSlideUniform()             );
        addType( new Move_SubtreeScale()                     );
        addType( new Move_SubtreeSwap()                      );
        addType( new Move_SPRNonclock()                      );
        addType( new Move_SpeciesNarrowExchange()            );
        addType( new Move_SpeciesNodeTimeSlideUniform()      );
        addType( new Move_SpeciesSubtreeScale()              );
        addType( new Move_SpeciesSubtreeScaleBeta()          );
        addType( new Move_SpeciesTreeScale()                 );
        addType( new Move_TreeScale()                        );
        addType( new Move_WeightedBranchLengthScale()        );

        /* Moves on character histories / data augmentation */
        addType( new Move_CharacterHistory() );
        // addType( new Move_NodeCharacterHistoryRejectionSample() );
        // addType( new Move_PathCharacterHistoryRejectionSample() );

        
        addType( new Move_IndependentPriorSampler<Real>( ) );
        addType( new Move_IndependentPriorSampler<RealPos>( ) );
        addType( new Move_IndependentPriorSampler<Natural>( ) );
        addType( new Move_IndependentPriorSampler<Integer>( ) );
        addType( new Move_IndependentPriorSampler<Probability>( ) );
        addType( new Move_IndependentPriorSampler<Simplex>( ) );
        addType( new Move_IndependentPriorSampler<ModelVector<RealPos> >( ) );
        addType( new Move_IndependentPriorSampler<ModelVector<Real> >( ) );
        addType( new Move_IndependentPriorSampler<ModelVector<Natural> >( ) );
//        addType( new Move_IndependentPriorSampler<RateGenerator>( ) );
        addType( new Move_IndependentPriorSampler<Tree>( ) );

        addType( new Move_ResampleFBD()                      );
    }
    catch(RbException& rbException)
    {

        RBOUT("Caught an exception while initializing moves in the workspace\n");
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
