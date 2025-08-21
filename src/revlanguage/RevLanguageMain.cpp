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
#include "RbSettings.h"

#ifdef RB_MPI
#include <mpi.h>
#endif

RevLanguageMain::RevLanguageMain(bool q)
    : quiet(q)
{

}


void RevLanguageMain::startRevLanguageEnvironment(const std::vector<std::string> &expressions, const std::optional<std::string>& filename, const std::vector<std::string> &args)
{
    auto& settings = RbSettings::userSettings();

    int rank = 0;
#ifdef RB_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    // 1. Load the modules
    try
    {
        RevLanguage::ModuleSystem::getModuleSystem().loadModules( settings.getModuleDir() );
    }    
    catch (RbException &e)
    {
        if (rank == 0) std::cout << e.getMessage() << std::endl;
    }


    // 2. Maybe show a header
    bool script_or_expr = filename or not expressions.empty();
    if (not script_or_expr and not quiet)
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
        int result = RevLanguage::Parser::getParser().processCommand(command_line, RevLanguage::Workspace::userWorkspacePtr());

        // We just hope for better input next time
        if ( result == 2 and not settings.getContinueOnError() )
        {
            RevClient::shutdown();
                
            exit(1);
        }
    }

    // 5. Evaluate any expressions given with -e expr
    for (auto expression: expressions)
    {
        // Should we be using RBOUT here?  It looks weird with the 2 spaces of padding.
        if (settings.getEcho() and rank == 0)
            std::cerr<<"> "<<expression<<"\n";

        int result = RevLanguage::Parser::getParser().processCommand(expression, RevLanguage::Workspace::userWorkspacePtr());
        
        // We just hope for better input next time
        if (result == 2 and not settings.getContinueOnError())
        {
            RevClient::shutdown();
                
            exit(1);
        }
    }

    // 6. Evaluate a filename if given.
    try
    {
        if (filename)
            RevClient::execute_file(*filename);
    }
    catch (const RbException& e)
    {
        // Try to give the same error messages as in RevLanguage::Parser::Execute( ) in revlanguage/parser/Parser.cpp
        if (rank == 0)
        {
            std::ostringstream msg;
            e.print(msg);
            RBOUT(msg.str());
        }
        std::exit(1);
    }
    catch (const std::exception& e)
    {
        if (rank == 0)
        {
            RBOUT(e.what());
        }
        std::exit(1);
    }
    catch (...)
    {
        if (rank == 0)
        {
            RBOUT("Error:\tunknown exception!");
        }
        std::exit(1);
    }
    
}

