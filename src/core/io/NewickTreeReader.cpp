#include <fstream>
#include <string>
#include <vector>

#include "NewickConverter.h"
#include "NewickTreeReader.h"
#include "RbException.h"
#include "RbFileManager.h"

namespace RevBayesCore { class Tree; }

using namespace RevBayesCore;


/**
 * Default constructor.
 * 
 * The default constructor does nothing except allocating the object.
 */
NewickTreeReader::NewickTreeReader() 
{
    
}



/**
 *
 */
std::vector<Tree*>* NewickTreeReader::readBranchLengthTrees(const path& fn)
{
    /* Open file */
    std::ifstream inFile( fn.string() );
    
    if ( !inFile )
        throw RbException()<<"Could not open file "<<fn;
    
    /* Read the whole file into a string */
    std::string input;
    bool firstline = true;
    while ( inFile.good() ) 
    {
        // Read a line
        std::string line;
        safeGetline( inFile, line );
        
        if (not firstline)
            input += "\n";
        else
            firstline = false;

        input += line;
    }
    
    return new std::vector<Tree*>(readNewicks(input));
}
