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
    
    /* Initialize */
    std::vector<Tree*>* trees = new std::vector<Tree*>();
    std::string commandLine;
    
    /* line-processing loop */
    while ( inFile.good() ) 
    {
        
        // Read a line
        std::string line;
        safeGetline( inFile, line );
        
        // skip empty lines
        if (line.length() == 0) 
        {
            continue;
        }
                
        
        NewickConverter c;
        Tree *blTree = c.convertFromNewick( line );
        
        trees->push_back( blTree );
    }
    
    return trees;
}
