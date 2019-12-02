#include "DelimitedCharacterDataReader.h"
#include "RbFileManager.h"
#include "RbConstants.h"

#include <boost/algorithm/string/trim_all.hpp>


using namespace RevBayesCore;


DelimitedCharacterDataReader::DelimitedCharacterDataReader(const std::string &fn, char d, size_t num_skipped) : DelimitedDataReader(fn, d, num_skipped)
{
    
    filename = fn;

    data = std::vector<std::vector<std::string> >( chars.size(), std::vector<std::string>() );
    for (size_t i = 0; i < chars.size(); ++i)
    {
        
        std::string name = chars[i][0];
        
        names.push_back( name );
        
        std::vector<std::string>& characters = data[i];
        for (size_t j = 1; j < chars[i].size(); ++j)
        {
            std::string s = chars[i][j];
            boost::trim( s );
            characters.push_back( s );
        }
        
    }
    
}


const std::vector<std::string>& DelimitedCharacterDataReader::getNames(void)
{
    return names;
}


const std::vector<std::vector<std::string> >& DelimitedCharacterDataReader::getData(void)
{
    return data;
}
