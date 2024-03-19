#include "RelativeNodeAgeWeightedConstraintsReader.h"

#include <cstdlib>
#include <string>

#include "RbException.h"


using namespace RevBayesCore;


RelativeNodeAgeWeightedConstraintsReader::RelativeNodeAgeWeightedConstraintsReader(const std::string &fn, std::string d, size_t ns, double thresh) : DelimitedDataReader(fn, d, ns)
{
    
    for (size_t i = 0; i < chars.size(); ++i)
    {
        
        if ( chars[i].size() < 5 )
        {
            throw RbException("Error while reading relative node age constraints. We expect 4 tip names and a weight per line and found fewer.");

        }
        if ( chars[i].size() > 5 )
        {
            throw RbException("Error while reading relative node age constraints. We expect 4 tip names and a weight per line and found more.");
            
        }
        //Get the two tip names defining the first (older) node:
        std::string name1 = chars[i][0];
        std::string name2 = chars[i][1];
        //Get the two tip names defining the second (younger) node:
        std::string name3 = chars[i][2];
        std::string name4 = chars[i][3];
        double weight = atof( chars[i][4].c_str() );
        if (weight >= thresh)
        {
            OlderYoungerConstraints.push_back(std::make_pair( std::make_pair( std::make_pair( name1, name2 ), std::make_pair( name3, name4 ) ), weight ) );
        }
    }
}


const std::vector< std::pair < std::pair < std::pair<std::string, std::string >, std::pair<std::string, std::string > >, double > >& RelativeNodeAgeWeightedConstraintsReader::getConstraints(void)
{
    return OlderYoungerConstraints;
}


const size_t RelativeNodeAgeWeightedConstraintsReader::getNumberOfConstraints(void) {
    
    return OlderYoungerConstraints.size();
}


