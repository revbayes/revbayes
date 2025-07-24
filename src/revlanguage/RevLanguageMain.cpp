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

RevLanguageMain::RevLanguageMain(bool b, bool b2) : batch_mode(b), show_header(b2)
{

}


void RevLanguageMain::startRevLanguageEnvironment(const std::vector<std::string> &expressions, const std::optional<std::string>& filename, const std::vector<std::string> &args)
{
    
    int pid = 0;
#ifdef RB_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
#endif
    
    // 1. Load the modules
    try
    {
        RevLanguage::ModuleSystem::getModuleSystem().loadModules( RbSettings::userSettings().getModuleDir() );
    }    
    catch (RbException &e)
    {
        std::cout << e.getMessage() << std::endl;
    }


    // 2. Maybe show a header
    if (show_header)
    {
        // Print a nifty message
        RbVersion version = RbVersion();
        RevLanguage::UserInterface::userInterface().output(version.getHeader(), false);
        RevLanguage::UserInterface::userInterface().output("", false);
    }
    
    // 3. Initialize the global workspace
    RevLanguage::Workspace::globalWorkspace().initializeGlobalWorkspace();

    // process the command line arguments as source file names    
    std::string line;
    std::string command_line;
    int result = 0;

    
    // 4. Set the args vector.

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

    // 5. Evaluate any expressions given with -e expr
    for (auto expression: expressions)
    {
        result = RevLanguage::Parser::getParser().processCommand(expression, RevLanguage::Workspace::userWorkspacePtr());
        
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

    // 6. Evaluate a filename if given.
    try
    {
        if (filename)
            RevClient::execute_file(*filename, false, true);
    }
    catch (const RbException& e)
    {
        RBOUT(e.getMessage());
    }
    catch (const std::exception& e)
    {
        if (rank == 0)
        {
            RBOUT(e.what());
        }
    }
    catch (...)
    {
        RBOUT("Error: unknown exception!");
    }
    
}

