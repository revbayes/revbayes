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


/// Monitors ///

/* Monitor types (in folder "monitors) */
#include "RlMonitor.h"
#include "Mntr_AncestralState.h"
#include "Mntr_JointConditionalAncestralState.h"
#include "Mntr_StochasticCharacterMapping.h"
#include "Mntr_File.h"
#include "Mntr_NexusFile.h"
#include "Mntr_ExtendedNewickFile.h"
#include "Mntr_Model.h"
#include "Mntr_Probability.h"
#include "Mntr_Screen.h"
#include "Mntr_StochasticVariable.h"
#include "Mntr_StochasticBranchRate.h"
#include "Mntr_StochasticBranchStateTimes.h"
#include "Mntr_CharacterHistoryNewickFile.h"
#include "Mntr_CharacterHistoryNhxFile.h"
#include "Mntr_CharacterHistorySummary.h"

/** Initialize global workspace */
void RevLanguage::Workspace::initializeMonitorGlobalWorkspace(void)
{
    
    try
    {
        ////////////////////////////////////////////////////////////////////////////////
        /* Add monitors (in folder "datatypes/inference/monitors") (alphabetic order) */
        ////////////////////////////////////////////////////////////////////////////////

		addTypeWithConstructor( new Mntr_AncestralState()                       );
        addTypeWithConstructor( new Mntr_JointConditionalAncestralState()       );
        addTypeWithConstructor( new Mntr_StochasticCharacterMapping()           );
        addTypeWithConstructor( new Mntr_ExtendedNewickFile()                   );
        addTypeWithConstructor( new Mntr_File()                                 );
        addTypeWithConstructor( new Mntr_NexusFile() );
        addTypeWithConstructor( new Mntr_Model()                                );
        addTypeWithConstructor( new Mntr_Probability()                          );
        addTypeWithConstructor( new Mntr_Screen()                               );
        addTypeWithConstructor( new Mntr_StochasticBranchRate()                 );
        addTypeWithConstructor( new Mntr_StochasticVariable()                   );
        addTypeWithConstructor( new Mntr_StochasticBranchStateTimes()           );
        addTypeWithConstructor( new Mntr_CharacterHistoryNewickFile()           );
        addTypeWithConstructor( new Mntr_CharacterHistoryNhxFile()              );
        addTypeWithConstructor( new Mntr_CharacterHistorySummary()              );

    }
    catch(RbException& rbException)
    {
        
        RBOUT("Caught an exception while initializing monitors in the workspace\n");
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


