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
#include <cmath>

/* Files including helper classes */
#include "RbException.h"
#include "RlUserInterface.h"
#include "Workspace.h"

/// Miscellaneous types ///

#include "DagMemberFunction.h"                        // for DagMemberFunction
#include "DeterministicNode.h"                        // for DeterministicNode
#include "DynamicNode.h"                              // for DynamicNode
#include "Func_workspaceVector.h"                     // for Func_workspaceV...
#include "RbVector.h"                                 // for RbVector
#include "RbVectorImpl.h"                             // for RbVectorImpl
#include "RevPtr.h"                                   // for RevPtr
#include "RlDagMemberFunction.h"                      // for DagMemberFunction
#include "RlDeterministicNode.h"                      // for DeterministicNode
#include "RlStochasticNode.h"                         // for StochasticNode
#include "RlStoppingRule.h"                           // for StoppingRule
#include "RlTypedDistribution.h"                      // for TypedDistribution
#include "RlTypedFunction.h"                          // for TypedFunction
#include "StochasticNode.h"                           // for StochasticNode
#include "Tree.h"                                     // for Tree, operator<<
#include "TypedDagNode.h"                             // for TypedDagNode
#include "TypedDistribution.h"                        // for TypedDistribution
#include "TypedFunction.h"                            // for TypedFunction
#include "WorkspaceToCoreWrapperObject.h"             // for WorkspaceToCore...

/* Base types (in folder "datatypes") */
#include "RevObject.h"

/* Container types (in folder "datatypes/container") */
#include "RlCorrespondenceAnalysis.h"

/* Container types (in folder "datatypes/math") */
#include "WorkspaceVector.h"

/* Container types (in folder "distributions/phylogenetics") */

/* Evolution types (in folder "datatypes/phylogenetics") */

/* Character state types (in folder "datatypes/phylogenetics/character") */

/* Character data types (in folder "datatypes/phylogenetics/datamatrix") */

/* Tree types (in folder "datatypes/phylogenetics/trees") */
#include "RlClade.h"

/* Taxon types (in folder "datatypes/phylogenetics") */
#include "RlTaxon.h"

/* Inference types (in folder "analysis") */
#include "RlBootstrapAnalysis.h"
#include "RlBurninEstimationConvergenceAssessment.h"
#include "RlHillClimber.h"
#include "RlMcmc.h"
#include "RlMcmcmc.h"
#include "RlModel.h"
#include "RlPathSampler.h"
#include "RlPosteriorPredictiveAnalysis.h"
#include "RlPosteriorPredictiveSimulation.h"
#include "RlPowerPosteriorAnalysis.h"
#include "RlSteppingStoneSampler.h"
#include "RlValidationAnalysis.h"
#include "RlTrace.h"
#include "RlTraceTree.h"
#include "RlAncestralStateTrace.h"

/// Stopping Rules ///
#include "RlMaxIterationStoppingRule.h"
#include "RlMaxTimeStoppingRule.h"
#include "RlMinEssStoppingRule.h"
#include "RlGelmanRubinStoppingRule.h"
#include "RlGewekeStoppingRule.h"
#include "RlStationarityStoppingRule.h"


/// Types ///

/* These types are needed as template types for the moves */
#include "RlMonitor.h"
#include "RlMove.h"
#include "RlTimeTree.h"

#include "Func_VectorMonitors.h"
#include "Func_VectorMoves.h"

#include "Probability.h"



/** Initialize global workspace */
void RevLanguage::Workspace::initializeTypeGlobalWorkspace(void)
{

    try
    {

        addFunction( new Func_VectorMoves()  );
        addFunction( new Func_VectorMonitors()  );

        addFunction( new Func_workspaceVector<Model>() );
        addFunction( new Func_workspaceVector<Trace>() );
        addFunction( new Func_workspaceVector<TraceTree>() );
		addFunction( new Func_workspaceVector<AncestralStateTrace>() );
        addFunction( new Func_workspaceVector<Monitor>() );
        addFunction( new Func_workspaceVector<Move>() );
        addFunction( new Func_workspaceVector<StoppingRule>() );

        /* Add math types (in folder "datatypes/math") */
        addType( new CorrespondenceAnalysis()                    );

        /* Add inference types (in folder "datatypes/inference") (alphabetic order) */
        addType( new BootstrapAnalysis()                             );
        addType( new BurninEstimationConvergenceAssessment()         );
        addType( new HillClimber()                                   );
        addType( new Mcmc()                                          );
        addType( new Mcmcmc()                                        );
        addType( new Model()                                         );
        addType( new PathSampler()                                   );
        addType( new PosteriorPredictiveAnalysis()                   );
        addType( new PosteriorPredictiveSimulation()                 );
        addType( new PowerPosteriorAnalysis()                        );
        addType( new SteppingStoneSampler()                          );
        addType( new ValidationAnalysis()                            );

        /* Add stopping rules (in folder "analysis/stoppingRules") (alphabetic order) */
        addType( new GelmanRubinStoppingRule()                   );
        addType( new GewekeStoppingRule()                        );
        addType( new MaxIterationStoppingRule()                  );
        addType( new MaxTimeStoppingRule()                       );
        addType( new MinEssStoppingRule()                        );
        addType( new StationarityStoppingRule()                  );

    }
    catch(RbException& rbException)
    {

        RBOUT("Caught an exception while initializing types in the workspace\n");
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
