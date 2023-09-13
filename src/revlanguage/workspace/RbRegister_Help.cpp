#include <cstdio>
#include <cstdlib>
#include <ostream>

#include "RbException.h"
#include "RbHelpSystem.h"
#include "RbHelpType.h"
#include "RlUserInterface.h"
#include "Workspace.h"

// The types of which we add extra help
#include "RlTimeTree.h"


/** 
 * Initialize extra help functions 
 */
void RevLanguage::Workspace::initializeExtraHelp(void)
{
    
    try
    {
        // add the help entry for this type to the global help system instance
//        MonteCarloAnalysis mca;
//        RevBayesCore::RbHelpSystem::getHelpSystem().addHelpType( static_cast<RevBayesCore::RbHelpType*>(mca.getHelpEntry()) );

        TimeTree tt;
        RevBayesCore::RbHelpSystem::getHelpSystem().addHelpType( static_cast<RevBayesCore::RbHelpType*>(tt.getHelpEntry()) );

    }
    catch(RbException& rbException)
    {
        
        RBOUT("Caught an exception while initializing the help system\n");
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


