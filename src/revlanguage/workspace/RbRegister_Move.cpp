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
#include <set>
#include <cstdlib>

/* Files including helper classes */
#include "AddWorkspaceVectorType.h"
#include "AddVectorizedWorkspaceType.h"
#include "RbException.h"
#include "RlUserInterface.h"
#include "Workspace.h"

/// Types ///

/* These types are needed as template types for the moves */
#include "RlBranchLengthTree.h"
#include "RlRateGenerator.h"


#include "Probability.h"
#include "RlSimplex.h"


/// Moves ///

/* Move types (in folder "datatypes/inference/moves") (grouped by parameter type) */
#include "RlMove.h"

/* Moves on real values */
#include "Move_Mirror.h"
#include "Move_MirrorMultiplier.h"
#include "Move_HSRFHyperpriorsGibbs.h"
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
#include "Move_EllipticalSliceSamplingLognormalIID.h"
#include "Move_EllipticalSliceSamplingSimple.h"
#include "Move_SynchronizedVectorFixedSingleElementSlide.h"
#include "Move_VectorBinarySwitch.h"
#include "Move_VectorSingleElementScale.h"
#include "Move_VectorSingleElementSlide.h"
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
#include "BetaSimplexProposal.h"


/* Moves on character histories/data augmentation */
// #include "Move_NodeCharacterHistoryRejectionSample.h"
// #include "Move_PathCharacterHistoryRejectionSample.h"
#include "Move_CharacterHistory.h"


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

/* Moves on continuous phyloprocesses (Brownian, multivariate Brownian, etc) */

/* Tree proposals (in folder "datatypes/inference/moves/tree") */
#include "Move_AddRemoveTip.h"
#include "Move_BurstEvent.h"
#include "Move_BranchLengthScale.h"
#include "Move_CollapseExpandFossilBranch.h"
#include "Move_IndependentTopology.h"
#include "Move_EmpiricalTree.h"
#include "Move_FNPR.h"
#include "Move_TipTimeSlideUniform.h"
#include "Move_GibbsPruneAndRegraft.h"
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
#include "Move_RootTimeScaleBactrian.h"
#include "Move_RootTimeSlideUniform.h"
#include "Move_SpeciesNarrowExchange.h"
#include "Move_SpeciesNodeTimeSlideUniform.h"
#include "Move_SpeciesSubtreeScale.h"
#include "Move_SpeciesSubtreeScaleBeta.h"
#include "Move_SpeciesTreeScale.h"
#include "Move_SubtreeScale.h"
#include "Move_SPRNonclock.h"
#include "Move_TreeScale.h"
#include "Move_WeightedNodeTimeSlide.h"


#include "Move_NarrowExchangeRateMatrix.h"

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
        addTypeWithConstructor( new Move_Mirror() );
        addTypeWithConstructor( new Move_MirrorMultiplier() );
        addTypeWithConstructor( new Move_Scale() );
        addTypeWithConstructor( new Move_GammaScale() );
        addTypeWithConstructor( new Move_RandomDive() );
        addTypeWithConstructor( new Move_ScaleBactrian() );
        addTypeWithConstructor( new Move_ScaleBactrianCauchy() );
        addTypeWithConstructor( new Move_Slide() );
        addTypeWithConstructor( new Move_SlideBactrian() );
        addTypeWithConstructor( new Move_SliceSampling() );
        addTypeWithConstructor( new Move_GMRFHyperpriorGibbs() );
        addTypeWithConstructor( new Move_GMRFUnevenGridHyperpriorGibbs() );
        /* Moves on probability */
        addTypeWithConstructor( new Move_BetaProbability() );

        /* compound moves */
//        addTypeWithConstructor("mvUpDownScale",         new Move_UpDownScale() );
        addTypeWithConstructor( new Move_AVMVN() );
        addTypeWithConstructor( new Move_UpDownTreeScale() );
        addTypeWithConstructor( new Move_UpDownSlide() );
        addTypeWithConstructor( new Move_UpDownSlideBactrian() );

		// compound moves on real values
        addTypeWithConstructor( new Move_LevyJumpSum() );
        addTypeWithConstructor( new Move_LevyJump() );

        /* Moves on integer values */
        addTypeWithConstructor( new Move_RandomIntegerWalk() );
        addTypeWithConstructor( new Move_RandomGeometricWalk() );
        addTypeWithConstructor( new Move_BinarySwitch() );

        /* Moves on simplices */
        addTypeWithConstructor( new Move_DirichletSimplex() );
        addTypeWithConstructor( new Move_BetaSimplex() );
        addTypeWithConstructor( new Move_ElementSwapSimplex() );

        /* Moves on vectors */
        addTypeWithConstructor( new Move_MultipleElementScale() );
        addTypeWithConstructor( new Move_HSRFHyperpriorsGibbs() );
        addTypeWithConstructor( new Move_HSRFUnevenGridHyperpriorsGibbs() );
        addTypeWithConstructor( new Move_SingleElementSlide() );
        addTypeWithConstructor( new Move_SingleElementScale() );
        addTypeWithConstructor( new Move_ShrinkExpand() );
        addTypeWithConstructor( new Move_ShrinkExpandScale() );
        addTypeWithConstructor( new Move_VectorBinarySwitch() );
        addTypeWithConstructor( new Move_VectorScale() );
        addTypeWithConstructor( new Move_VectorSlide() );
        addTypeWithConstructor( new Move_VectorSlideRecenter() );
        addTypeWithConstructor( new Move_ElementScale() );
        addTypeWithConstructor( new Move_ElementSlide() );
        addTypeWithConstructor( new Move_VectorSingleElementScale() );
        addTypeWithConstructor( new Move_VectorSingleElementSlide() );
        addTypeWithConstructor( new Move_VectorFixedSingleElementSlide() );
        addTypeWithConstructor( new Move_EllipticalSliceSamplingLognormalIID() );
        addTypeWithConstructor( new Move_EllipticalSliceSamplingSimple() );
        addTypeWithConstructor( new Move_SynchronizedVectorFixedSingleElementSlide() );

        /* Moves on matrices of real values */
        addTypeWithConstructor( new Move_MatrixSingleElementScale() );
        addTypeWithConstructor( new Move_MatrixSingleElementSlide() );

        /* Moves on matrices of correlations */
        addTypeWithConstructor( new Move_CorrelationMatrixUpdate()                   );
        addTypeWithConstructor( new Move_CorrelationMatrixRandomWalk()               );
        addTypeWithConstructor( new Move_CorrelationMatrixSingleElementBeta()        );
        addTypeWithConstructor( new Move_CorrelationMatrixSpecificElementBeta()      );
        addTypeWithConstructor( new Move_CorrelationMatrixElementSwap()              );

//        addTypeWithConstructor( new Move_CorrelationMatrixPartialSingleElementBeta() );
//        addTypeWithConstructor( new Move_CorrelationMatrixReparameterization()       );
//        addTypeWithConstructor( new Move_CorrelationMatrixExpansion()                );

        /* Moves on matrices of real values */
        addTypeWithConstructor( new Move_MatrixRealSymmetricSlide() );

        /* Moves on matrices of real values */
        addTypeWithConstructor( new Move_ConjugateInverseWishart() );

        /* Moves of random adjacency graphs */
        addTypeWithConstructor( new Move_GraphFlipEdge() );
        addTypeWithConstructor( new Move_GraphFlipClique() );
        addTypeWithConstructor( new Move_GraphShiftEdge() );

        /* Moves on continuous character data (matrices of real values) */
        addTypeWithConstructor( new Move_ContinuousCharacterDataSlide() );

        /* Moves on mixtures (in folder "datatypes/inference/moves/mixture") */
        addTypeWithConstructor( new Move_DPPTableValueUpdate<RealPos>( new RevBayesCore::ScaleProposal( NULL, 1.0 ) ) );
        addTypeWithConstructor( new Move_DPPTableValueUpdate<Simplex>( new RevBayesCore::BetaSimplexProposal( NULL, 10.0 ) ) );

//        addTypeWithConstructor("mvDPPScaleCatVals",                new Move_DPPScaleCatValsMove() );
//        addTypeWithConstructor("mvDPPScaleCatAllocateAux",         new Move_DPPScaleCatAllocateAux() );
        addTypeWithConstructor( new Move_DPPAllocateAuxGibbsMove<Real>() );
        addTypeWithConstructor( new Move_DPPAllocateAuxGibbsMove<RealPos>() );
        addTypeWithConstructor( new Move_DPPAllocateAuxGibbsMove<Probability>() );
        addTypeWithConstructor( new Move_DPPAllocateAuxGibbsMove<Integer>() );
        addTypeWithConstructor( new Move_DPPAllocateAuxGibbsMove<Natural>() );
        addTypeWithConstructor( new Move_DPPAllocateAuxGibbsMove<Simplex>() );
        addTypeWithConstructor( new Move_DPPGibbsConcentration( ) );
        addTypeWithConstructor( new Move_MixtureAllocation<Real>( ) );
        addTypeWithConstructor( new Move_MixtureAllocation<RealPos>( ) );
        addTypeWithConstructor( new Move_MixtureAllocation<Natural>( ) );
        addTypeWithConstructor( new Move_MixtureAllocation<Integer>( ) );
        addTypeWithConstructor( new Move_MixtureAllocation<Probability>( ) );
        addTypeWithConstructor( new Move_MixtureAllocation<Simplex>( ) );
        addTypeWithConstructor( new Move_MixtureAllocation<RateGenerator>( ) );
        addTypeWithConstructor( new Move_MixtureAllocation<Tree>( ) );
        addTypeWithConstructor( new Move_GibbsMixtureAllocation<Real>( ) );
        addTypeWithConstructor( new Move_GibbsMixtureAllocation<RealPos>( ) );
        addTypeWithConstructor( new Move_GibbsMixtureAllocation<Natural>( ) );
        addTypeWithConstructor( new Move_GibbsMixtureAllocation<Integer>( ) );
        addTypeWithConstructor( new Move_GibbsMixtureAllocation<Probability>( ) );
        addTypeWithConstructor( new Move_GibbsMixtureAllocation<Simplex>( ) );
        addTypeWithConstructor( new Move_GibbsMixtureAllocation<RateGenerator>( ) );
        addTypeWithConstructor( new Move_UPPAllocation<RealPos>() );

        addTypeWithConstructor( new Move_ReversibleJumpSwitch<Real>( )                  );
        addTypeWithConstructor( new Move_ReversibleJumpSwitch<RealPos>( )               );
        addTypeWithConstructor( new Move_ReversibleJumpSwitch<Natural>( )               );
        addTypeWithConstructor( new Move_ReversibleJumpSwitch<Integer>( )               );
        addTypeWithConstructor( new Move_ReversibleJumpSwitch<Probability>( )           );
        addTypeWithConstructor( new Move_ReversibleJumpSwitch<Simplex>( )               );
        addTypeWithConstructor( new Move_ReversibleJumpSwitch<ModelVector<Natural> >( ) );
        addTypeWithConstructor( new Move_ReversibleJumpSwitch<Tree>( )                  );
        
        addTypeWithConstructor( new Move_MultiValueEventBirthDeath()                    );
        addTypeWithConstructor( new Move_MultiValueEventScale()                         );
        addTypeWithConstructor( new Move_MultiValueEventSlide()                         );



        addTypeWithConstructor( new Move_BirthDeathEventContinuous()                    );
        addTypeWithConstructor( new Move_BirthDeathEventDiscrete()                      );
        addTypeWithConstructor( new Move_ContinuousEventScale()                         );
        addTypeWithConstructor( new Move_DiscreteEventCategoryRandomWalk()              );
        addTypeWithConstructor( new Move_EventTimeBeta()                                );
        addTypeWithConstructor( new Move_EventTimeSlide()                               );
        addTypeWithConstructor( new Move_BirthDeathFromAgeEvent()                       );
        addTypeWithConstructor( new Move_GibbsDrawCharacterHistory()                    );

        /* Tree proposals (in folder "datatypes/inference/moves/tree") */
        addTypeWithConstructor( new Move_AddRemoveTip()                     );
        addTypeWithConstructor( new Move_BurstEvent()                       );
        addTypeWithConstructor( new Move_BranchLengthScale()                );
        addTypeWithConstructor( new Move_CollapseExpandFossilBranch()       );
        addTypeWithConstructor( new Move_IndependentTopology()              );
		addTypeWithConstructor( new Move_EmpiricalTree()                    );
        addTypeWithConstructor( new Move_FNPR()                             );
        addTypeWithConstructor( new Move_GibbsPruneAndRegraft()             );
        addTypeWithConstructor( new Move_LayeredScaleProposal()             );
        addTypeWithConstructor( new Move_NarrowExchange()                   );
        addTypeWithConstructor( new Move_NNIClock()                         );
        addTypeWithConstructor( new Move_NNINonclock()                      );
        addTypeWithConstructor( new Move_NodeRateTimeSlideUniform()         );
        addTypeWithConstructor( new Move_NodeTimeScale()                    );
        addTypeWithConstructor( new Move_NodeTimeSlidePathTruncatedNormal() );
        addTypeWithConstructor( new Move_NodeTimeSlideUniform()             );
        addTypeWithConstructor( new Move_NodeTimeSlideUniformAgeConstrained());
        addTypeWithConstructor( new Move_NodeTimeSlideBeta()                );
        addTypeWithConstructor( new Move_RateAgeBetaShift()                 );
        addTypeWithConstructor( new Move_RootTimeScaleBactrian()            );
        addTypeWithConstructor( new Move_RootTimeSlideUniform()             );
        addTypeWithConstructor( new Move_SubtreeScale()                     );
        addTypeWithConstructor( new Move_SPRNonclock()                      );
        addTypeWithConstructor( new Move_SpeciesNarrowExchange()            );
        addTypeWithConstructor( new Move_SpeciesNodeTimeSlideUniform()      );
        addTypeWithConstructor( new Move_SpeciesSubtreeScale()              );
        addTypeWithConstructor( new Move_SpeciesSubtreeScaleBeta()          );
        addTypeWithConstructor( new Move_TipTimeSlideUniform()              );
        addTypeWithConstructor( new Move_SpeciesTreeScale()                 );
        addTypeWithConstructor( new Move_TreeScale()                        );
        addTypeWithConstructor( new Move_NarrowExchangeRateMatrix()         );

        /* Moves on character histories / data augmentation */
        addTypeWithConstructor( new Move_CharacterHistory() );
        // addTypeWithConstructor( new Move_NodeCharacterHistoryRejectionSample() );
        // addTypeWithConstructor( new Move_PathCharacterHistoryRejectionSample() );

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
