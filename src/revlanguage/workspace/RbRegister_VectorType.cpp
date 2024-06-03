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
#include "AddWorkspaceVectorType.h"
#include "RlUserInterface.h"

/// Miscellaneous types ///

#include "AbstractHomologousDiscreteCharacterData.h"
#include "AbstractDiscreteTaxonData.h"
#include "AverageDistanceMatrix.h"
#include "CharacterHistoryRateModifier.h"
#include "Clade.h"
#include "CladogeneticProbabilityMatrix.h"
#include "CladogeneticSpeciationRateMatrix.h"
#include "ConstantNode.h"
#include "ContinuousCharacterData.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DistanceMatrix.h"
#include "DynamicNode.h"
#include "Func__vectorIndexOperator.h"
#include "Func_modelVector.h"
#include "Func_workspaceVector.h"
#include "IndirectReferenceFunction.h"
#include "MatrixBoolean.h"
#include "MatrixReal.h"
#include "ModelObject.h"
#include "RateGenerator.h"
#include "RbException.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "RevPtr.h"
#include "RlConstantNode.h"
#include "RlDeterministicNode.h"
#include "RlDistribution.h"
#include "RlSiteMixtureModel.h"
#include "RlTypedDistribution.h"
#include "RlTree.h"
#include "RlTypedFunction.h"
#include "Taxon.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"
#include "VectorFunction.h"
#include "VectorIndexOperator.h"
#include "Workspace.h"
#include "WorkspaceToCoreWrapperObject.h"

/* Base types (in folder "datatypes") */
#include "RevObject.h"

/* Container types (in folder "datatypes/container") */
#include "RlMatrixBoolean.h"
#include "RlMatrixReal.h"
#include "RlMatrixRealPos.h"
#include "RlMatrixRealSymmetric.h"
#include "RlStochasticMatrix.h"

/* Container types (in folder "datatypes/math") */
#include "ModelVector.h"
#include "WorkspaceVector.h"

/* Character data types (in folder "datatypes/phylogenetics/datamatrix") */
#include "RlAbstractDiscreteTaxonData.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlContinuousCharacterData.h"

/* Tree types (in folder "datatypes/phylogenetics/trees") */
#include "RlClade.h"

/* Taxon types (in folder "datatypes/phylogenetics") */
#include "RlCladogeneticProbabilityMatrix.h"
#include "RlCladogeneticSpeciationRateMatrix.h"
#include "RlAverageDistanceMatrix.h"
#include "RlDistanceMatrix.h"
#include "RlCharacterHistoryRateModifier.h"
#include "RlTaxon.h"

/* Regional feature types (in folder "datatypes/phylogenetics/geography") */
#include "RlRegionalFeatures.h"
#include "RlRegionalFeatureLayer.h"

/// Types ///

#include "Probability.h"

/* These types are needed as template types for the moves */
#include "RlBranchLengthTree.h"
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
        AddWorkspaceVectorType<AverageDistanceMatrix,3>::addTypeToWorkspace( *this, new AverageDistanceMatrix() );
        AddWorkspaceVectorType<DistanceMatrix,3>::addTypeToWorkspace( *this, new DistanceMatrix() );
        AddWorkspaceVectorType<MatrixBoolean,3>::addTypeToWorkspace( *this, new MatrixBoolean() );
        AddWorkspaceVectorType<MatrixReal,3>::addTypeToWorkspace( *this, new MatrixReal() );
        AddWorkspaceVectorType<MatrixRealPos,3>::addTypeToWorkspace( *this, new MatrixRealPos() );
        AddWorkspaceVectorType<MatrixRealSymmetric,3>::addTypeToWorkspace( *this, new MatrixRealSymmetric() );
        AddWorkspaceVectorType<StochasticMatrix,3>::addTypeToWorkspace( *this, new StochasticMatrix() );
        AddWorkspaceVectorType<AbstractDiscreteTaxonData,3>::addTypeToWorkspace( *this, new AbstractDiscreteTaxonData() );
        AddWorkspaceVectorType<AbstractHomologousDiscreteCharacterData,3>::addTypeToWorkspace( *this, new AbstractHomologousDiscreteCharacterData() );
        AddWorkspaceVectorType<ContinuousCharacterData,3>::addTypeToWorkspace( *this, new ContinuousCharacterData() );
        AddWorkspaceVectorType<CharacterHistoryRateModifier,3>::addTypeToWorkspace( *this, new CharacterHistoryRateModifier() );
        AddWorkspaceVectorType<TimeTree,3>::addTypeToWorkspace( *this, new TimeTree() );
        AddWorkspaceVectorType<BranchLengthTree,3>::addTypeToWorkspace( *this, new BranchLengthTree() );
        AddWorkspaceVectorType<Tree,3>::addTypeToWorkspace( *this, new Tree() );
        AddWorkspaceVectorType<Clade,3>::addTypeToWorkspace( *this, new Clade() );
        AddWorkspaceVectorType<SiteMixtureModel,3>::addTypeToWorkspace( *this, new SiteMixtureModel() );
//        AddWorkspaceVectorType<RlRegionalFeatures,3>::addTypeToWorkspace( *this, new RlRegionalFeatures() );
        AddWorkspaceVectorType<RlRegionalFeatureLayer,3>::addTypeToWorkspace( *this, new RlRegionalFeatureLayer() );
        //        AddWorkspaceVectorType<Dist_bdp,3>::addTypeToWorkspace( *this, new Dist_bdp() );

//        AddWorkspaceVectorType<Dist_unif,1>::addTypeToWorkspace( *this, new Dist_unif() );
//        this->addFunction(new Func_workspaceVector<Dist_unif>() );
        addFunction(new Func_workspaceVector<Distribution>() );
        addFunction( new Func_workspaceVector<TypedDistribution<Natural> >() );
        addFunction(new Func_workspaceVector<TypedDistribution<Real> >() );
        addFunction(new Func_workspaceVector<TypedDistribution<RealPos> >() );
        addFunction(new Func_workspaceVector<TypedDistribution<Probability> >() );
        addFunction(new Func_workspaceVector<TypedDistribution<ModelVector<Natural> > >() );
        addFunction(new Func_workspaceVector<TypedDistribution<ModelVector<Real> > >() );
        addFunction(new Func_workspaceVector<TypedDistribution<ModelVector<RealPos> > >() );
        addFunction(new Func_workspaceVector<TypedDistribution<ModelVector<Probability> > >() );
        addFunction( new Func_workspaceVector<TypedDistribution<TimeTree> >() );
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
