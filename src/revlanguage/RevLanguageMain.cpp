#include "ModuleSystem.h"
#include "RevLanguageMain.h"
#include "Parser.h"
#include "RbException.h"
#include "RbFileManager.h"
#include "RbSettings.h"
#include "Workspace.h"
#include "RlUserInterface.h"
#include "RbVersion.h"
#include <iostream>
#include <string>
#include <cstdlib>

#ifdef RB_MPI
#include <mpi.h>
#endif

RevLanguageMain::RevLanguageMain(bool b) : batch_mode(b)
{

}


void RevLanguageMain::startRevLanguageEnvironment(std::vector<std::string> sourceFiles)
{
    
    int pid = 0;
#ifdef RB_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
#endif
    
    // load the modules
    try
    {
        RevLanguage::ModuleSystem::getModuleSystem().loadModules( RbSettings::userSettings().getModuleDir() );
    }    
    catch (RbException &e)
    {
        std::cout << e.getMessage() << std::endl;
    }


    // Print a nifty message
    RbVersion version = RbVersion();
    RevLanguage::UserInterface::userInterface().output(version.getHeader(), false);
    RevLanguage::UserInterface::userInterface().output("", false);
    
    RevLanguage::Workspace::globalWorkspace().initializeGlobalWorkspace();

    // process the command line arguments as source file names    
    std::string line;
    std::string commandLine;
    int result = 0;

    for (unsigned int i =0 ; i < sourceFiles.size(); i++)
    {
        line = "source(\"" + sourceFiles[i] + "\")";
        
        // let only the master process print to the screen
        if ( pid == 0 )
        {
            std::cout << "> " << line << std::endl;
        }
        
        // Process the command line
        if (result == 1)
        {
            commandLine += line;
        }
        else
        {
            commandLine = line;
        }
        
        result = RevLanguage::Parser::getParser().processCommand(commandLine, &RevLanguage::Workspace::userWorkspace());

        // We just hope for better input next time
        if (result == 2)
        {
            result = 0;

            if( batch_mode == true )
            {
                RevLanguage::Workspace::userWorkspace().clear();
                RevLanguage::Workspace::globalWorkspace().clear();
                
#ifdef RB_MPI
                MPI_Finalize();
#endif

                exit(1);
            }
        }
        
    }
    
}

