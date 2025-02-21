#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>

#include "ModuleSystem.h"
#include "RevLanguageMain.h"
#include "Parser.h"
#include "RbException.h"
#include "RbSettings.h"
#include "Workspace.h"
#include "RlUserInterface.h"
#include "RbVersion.h"
#include "StringUtilities.h"
#include "RevClient.h"        // for RevClient::shutdown()

#ifdef RB_MPI
#include <mpi.h>
#endif

RevLanguageMain::RevLanguageMain(bool b) : batch_mode(b)
{

}


void RevLanguageMain::startRevLanguageEnvironment(const std::vector<std::string> &args, const std::vector<std::string> &source_files)
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
    std::string command_line;
    int result = 0;

    // Ensure that args exists and has size 0.
    command_line = "args = [\"\"]; args.erase(\"\")";
    RevLanguage::Parser::getParser().processCommand(command_line, RevLanguage::Workspace::userWorkspacePtr());

    for (unsigned int i =0 ; i < args.size(); ++i)
    {
        if ( StringUtilities::isNumber( args[i] ) )
        {
            command_line = "args[" + StringUtilities::to_string(i+1) + "] = " + args[i];
        }
        else
        {
            command_line = "args[" + StringUtilities::to_string(i+1) + "] = \"" + args[i] + "\"";
        }
        result = RevLanguage::Parser::getParser().processCommand(command_line, RevLanguage::Workspace::userWorkspacePtr());
        
        // We just hope for better input next time
        if (result == 2)
        {
            result = 0;
            
            if( batch_mode == true )
            {
                RevClient::shutdown();
                
                exit(1);
            }
        }
    }
    
    for (unsigned int i =0 ; i < source_files.size(); ++i)
    {
        line = "source(\"" + source_files[i] + "\")";
        
        // let only the master process print to the screen
        if ( pid == 0 )
        {
            std::cout << "> " << line << std::endl;
        }
        
        // Process the command line
        if (result == 1)
        {
            command_line += line;
        }
        else
        {
            command_line = line;
        }
        
        result = RevLanguage::Parser::getParser().processCommand(command_line, RevLanguage::Workspace::userWorkspacePtr());

        // We just hope for better input next time
        if (result == 2)
        {
            result = 0;

            if( batch_mode == true )
            {
                RevClient::shutdown();

                exit(1);
            }
        }
        
    }
    
}

