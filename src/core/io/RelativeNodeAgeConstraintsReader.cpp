#include "RelativeNodeAgeConstraintsReader.h"
#include "RbException.h"
#include "RbFileManager.h"

#include "RbConstants.h"


using namespace RevBayesCore;


RelativeNodeAgeConstraintsReader::RelativeNodeAgeConstraintsReader(const std::string &fn, char d, size_t ns) : DelimitedDataReader(fn, d, ns)
{
    
    filename = fn;
    
    for (size_t i = 0; i < chars.size(); ++i)
    {
        
        if ( chars[i].size() < 4 )
        {
            throw RbException("Error while reading relative node age constraints. We expect 4 tip names per line and found fewer.");

        }
        if ( chars[i].size() > 4 )
        {
            throw RbException("Error while reading relative node age constraints. We expect 4 tip names per line and found more.");
            
        }
        //Get the two tip names defining the first (older) node:
        std::string name1 = chars[i][0];
        std::string name2 = chars[i][1];
        //Get the two tip names defining the second (younger) node:
        std::string name3 = chars[i][2];
        std::string name4 = chars[i][3];
        OlderYoungerConstraints.push_back(std::make_pair( std::make_pair( name1, name2 ), std::make_pair( name3, name4 ) ) );
        
    }
}


const std::vector< std::pair < std::pair<std::string, std::string >, std::pair<std::string, std::string > > >& RelativeNodeAgeConstraintsReader::getConstraints(void)
{
    return OlderYoungerConstraints;
}


const size_t RelativeNodeAgeConstraintsReader::getNumberOfConstraints(void) {
    
    return OlderYoungerConstraints.size();
}


