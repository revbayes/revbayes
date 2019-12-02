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
#include "RlUserInterface.h"

/// Miscellaneous types ///

/* Base types (in folder "datatypes") */
#include "RevObject.h"
#include "AbstractModelObject.h"

/* Container types (in folder "datatypes/container") */
#include "RlCorrespondenceAnalysis.h"
#include "RlMatrixReal.h"
#include "RlMatrixRealPos.h"
#include "RlMatrixRealSymmetric.h"
#include "RlRateGeneratorSequence.h"
#include "RlRateMatrix.h"
#include "RlSimplex.h"

/* Container types (in folder "datatypes/math") */
#include "ModelVector.h"
#include "WorkspaceVector.h"

/* Character data types (in folder "datatypes/phylogenetics/datamatrix") */
#include "RlAbstractCharacterData.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlContinuousCharacterData.h"

/* Tree types (in folder "datatypes/phylogenetics/trees") */
#include "RlClade.h"
#include "RlRootedTripletDistribution.h"


/* Taxon types (in folder "datatypes/phylogenetics") */
#include "RlCladogeneticProbabilityMatrix.h"
#include "RlCladogeneticSpeciationRateMatrix.h"
#include "RlDistanceMatrix.h"
#include "RlCharacterHistoryRateModifier.h"
#include "RlTaxon.h"

#include "Dist_unif.h"

/// Types ///

/* These types are needed as template types for the moves */
#include "RlBranchLengthTree.h"
#include "RlMonitor.h"
#include "RlMove.h"
#include "RlRateGenerator.h"
#include "RlTimeTree.h"

/** Initialize global workspace */
void RevLanguage::Workspace::initializeVectorTypeGlobalWorkspace(void)
{
    
    try
    {
        
        AddWorkspaceVectorType<Taxon,4>::addTypeToWorkspace( *this, new Taxon() );
        AddWorkspaceVectorType<RateGenerator,3>::addTypeToWorkspace( *this, new RateGenerator() );
        AddWorkspaceVectorType<CladogeneticProbabilityMatrix,3>::addTypeToWorkspace( *this, new CladogeneticProbabilityMatrix() );
        AddWorkspaceVectorType<CladogeneticSpeciationRateMatrix,3>::addTypeToWorkspace( *this, new CladogeneticSpeciationRateMatrix() );
        AddWorkspaceVectorType<DistanceMatrix,3>::addTypeToWorkspace( *this, new DistanceMatrix() );
        AddWorkspaceVectorType<MatrixReal,3>::addTypeToWorkspace( *this, new MatrixReal() );
        AddWorkspaceVectorType<MatrixRealPos,3>::addTypeToWorkspace( *this, new MatrixRealPos() );
        AddWorkspaceVectorType<MatrixRealSymmetric,3>::addTypeToWorkspace( *this, new MatrixRealSymmetric() );
        AddWorkspaceVectorType<AbstractHomologousDiscreteCharacterData,3>::addTypeToWorkspace( *this, new AbstractHomologousDiscreteCharacterData() );
        AddWorkspaceVectorType<ContinuousCharacterData,3>::addTypeToWorkspace( *this, new ContinuousCharacterData() );
        AddWorkspaceVectorType<CharacterHistoryRateModifier,3>::addTypeToWorkspace( *this, new CharacterHistoryRateModifier() );
        AddWorkspaceVectorType<TimeTree,3>::addTypeToWorkspace( *this, new TimeTree() );
        AddWorkspaceVectorType<BranchLengthTree,3>::addTypeToWorkspace( *this, new BranchLengthTree() );
        AddWorkspaceVectorType<Tree,3>::addTypeToWorkspace( *this, new Tree() );
        AddWorkspaceVectorType<Clade,3>::addTypeToWorkspace( *this, new Clade() );
        //        AddWorkspaceVectorType<Dist_bdp,3>::addTypeToWorkspace( *this, new Dist_bdp() );

//        AddWorkspaceVectorType<Dist_unif,1>::addTypeToWorkspace( *this, new Dist_unif() );
//        this->addFunction(new Func_workspaceVector<Dist_unif>() );
        this->addFunction(new Func_workspaceVector<Distribution>() );
//        this->addFunction(new Func_workspaceVector<ContinuousDistribution>() );
//        AddWorkspaceVectorType<Distribution,1>::addTypeToWorkspace( *this, new Distribution() );

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
