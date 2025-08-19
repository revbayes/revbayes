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
#include <string>

/* Files including helper classes */
#include "AddWorkspaceVectorType.h"
#include "RbException.h"
#include "RevAbstractType.h"
#include "RlUserInterface.h"
#include "Workspace.h"

/// Miscellaneous types ///

#include "ConstantNode.h"               // for ConstantNode
#include "DagNode.h"                    // for DagNode
#include "DeterministicNode.h"          // for DeterministicNode
#include "DynamicNode.h"                // for DynamicNode
#include "Func__vectorIndexOperator.h"  // for Func__vectorIndexOperator
#include "Func_modelVector.h"           // for Func_modelVector
#include "IndirectReferenceFunction.h"  // for IndirectReferenceFunction
#include "ModelObject.h"                // for ModelObject, ModelObject<>::v...
#include "RbBoolean.h"                  // for Boolean, operator<<
#include "RbVector.h"                   // for RbVector
#include "RbVectorImpl.h"               // for RbVectorImpl
#include "RevNullObject.h"              // for RevNullObject
#include "RevPtr.h"                     // for RevPtr
#include "RlConstantNode.h"             // for ConstantNode
#include "RlDeterministicNode.h"        // for DeterministicNode
#include "RlTypedFunction.h"            // for TypedFunction
#include "Simplex.h"                    // for Simplex
#include "TypedDagNode.h"               // for TypedDagNode
#include "TypedFunction.h"              // for TypedFunction
#include "UserFunctionNode.h"           // for UserFunctionNode
#include "VectorFunction.h"             // for VectorFunction
#include "VectorIndexOperator.h"        // for VectorIndexOperator

/* Base types (in folder "datatypes") */
#include "RevObject.h"

/* Primitive types (in folder "datatypes/basic") */
#include "Integer.h"
#include "IntegerPos.h"
#include "Natural.h"
#include "Probability.h"
#include "RlBoolean.h"
#include "RlString.h"
#include "Real.h"
#include "RealPos.h"

/* Container types (in folder "datatypes/container") */

#include "RlSimplex.h"

/* Container types (in folder "datatypes/math") */
#include "ModelVector.h"


/** Initialize global workspace */
void RevLanguage::Workspace::initializeBasicTypeGlobalWorkspace(void)
{
    
    try
    {
        /* Add types: add a dummy variable which we use for type checking, conversion checking and other tasks. */
        
        /* Add base types (in folder "datatypes") */
        addType( new RevAbstractType( RevObject::getClassTypeSpec(), &RevNullObject::getInstance() ) );
        
        /* Add primitive types (in folder "datatypes/basic") (alphabetic order) */
        AddWorkspaceVectorType<Integer,4>::addTypeToWorkspace(     *this, new Integer()     );
        AddWorkspaceVectorType<Natural,4>::addTypeToWorkspace(     *this, new Natural()     );
        AddWorkspaceVectorType<IntegerPos,4>::addTypeToWorkspace(  *this, new IntegerPos()  );
        AddWorkspaceVectorType<Probability,4>::addTypeToWorkspace( *this, new Probability() );
        AddWorkspaceVectorType<Real,4>::addTypeToWorkspace(        *this, new Real()        );
        AddWorkspaceVectorType<RealPos,4>::addTypeToWorkspace(     *this, new RealPos()     );
        AddWorkspaceVectorType<RlBoolean,4>::addTypeToWorkspace(   *this, new RlBoolean()   );
        AddWorkspaceVectorType<RlString,4>::addTypeToWorkspace(    *this, new RlString()    );
        AddWorkspaceVectorType<Simplex,4>::addTypeToWorkspace(     *this, new Simplex()     );


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


