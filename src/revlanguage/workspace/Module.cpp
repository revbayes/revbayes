#include "Module.h"
#include "RbFileManager.h"

using namespace RevLanguage;


Module::Module(const std::string &fn)
{
    // first we need to load the file
    std::ifstream readStream(fn);
    
    /* File-processing loop */
    while ( readStream.good() )
    {
        // Read a line
        std::string line;
        RevBayesCore::safeGetline(readStream, line);
        
        commandLines.push_back( line );
    }
}


Module::~Module( void )
{
    
}


const std::vector<std::string>& Module::getCommandLines( void ) const
{
    return commandLines;
}

