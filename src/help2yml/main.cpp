#include <vector>
#include <iostream>
#include <fstream>

#include "RbFileManager.h"
#include "RbHelpDistribution.h"
#include "RbHelpMonitor.h"
#include "RbHelpMove.h"
#include "RbHelpSystem.h"
#include "StringUtilities.h"
#include "YAMLHelpRenderer.h"
#include "Workspace.h"

using namespace RevBayesCore;

int main(int argc, const char * argv[])
{
    path file = "help.yml";

    if( argc > 1 )
    {
        path prefix = std::string(argv[1]);

        if (not is_directory(prefix))
        {
            std::cerr << "Error: Directory " << prefix << " does not exist" << std::endl;
            exit(1);
        }
        file = prefix / file;
    }

    std::string function_entry_result, type_entry_result, dist_entry_result, mntr_entry_result, move_entry_result, tmp;
    
    // first we need to initialize the workspace
    RevLanguage::Workspace::globalWorkspace().initializeGlobalWorkspace();
    
    RevBayesCore::RbHelpSystem& help = RevBayesCore::RbHelpSystem::getHelpSystem();
    YAML::HelpRenderer renderer = YAML::HelpRenderer();


    std::cout << std::endl << "Generating help files..." << std::endl;

    std::fstream fs;
    fs.open(file.string(), std::fstream::out | std::fstream::trunc );

    const std::set<std::string> &functionEntryNames = help.getFunctionEntries();
    const std::set<std::string> &typeEntryNames = help.getTypeEntries();
    
    for (std::set<std::string>::const_iterator it=functionEntryNames.begin(); it!=functionEntryNames.end(); ++it)
    {
        std::string name = *it;
        
        std::string base = name;
        base.erase(base.begin());
        std::string base_upper = StringUtilities::firstCharToUpper(base);

        if ( typeEntryNames.find("dn"+base) != typeEntryNames.end() || typeEntryNames.find("dn"+base_upper) != typeEntryNames.end() )
        {
            continue;
        }

        const RevBayesCore::RbHelpFunction& entry = static_cast<const RevBayesCore::RbHelpFunction&>( help.getHelp( name ) );
        
        if (name.size() > 0 && entry.getUsage() != "c_name()")
        {
            std::cout << "Adding function with name:\t" << name << std::endl;

            fs << renderer.renderHelp(entry);
        }
        
    }
    
    for (std::set<std::string>::const_iterator it=typeEntryNames.begin(); it!=typeEntryNames.end(); ++it)
    {
        std::string name = *it;
        
        StringUtilities::replaceSubstring(name,"[]","");

        // skip vector types
        if( *it != name )
        {
            continue;
        }

        if (name.size() > 0 && name != "c_name" && functionEntryNames.find(name) == functionEntryNames.end() )
        {
            std::cerr << "Adding type with name:\t" << name << std::endl;
            
            const RevBayesCore::RbHelpType& typeEntry = static_cast<const RevBayesCore::RbHelpType&>( help.getHelp( name ) );

            fs << renderer.renderHelp(typeEntry);
            
            
        }
        
    }
    
    fs.close();

    std::cout << "\nHelp files are updated." << std::endl;

    return 0;
}
