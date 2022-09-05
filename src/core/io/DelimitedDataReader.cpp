#include "DelimitedDataReader.h"

#include <functional>
#include <algorithm>
#include <string>
#include <sstream> // IWYU pragma: keep

#include "RbFileManager.h"
#include "RbException.h"
#include "StringUtilities.h"

using namespace RevBayesCore;

DelimitedDataReader::DelimitedDataReader(const path &fn, std::string d, size_t lines_skipped) :
    filename(fn), 
    delimiter(d),
    chars()
{

    readData( lines_skipped );
    
}

void DelimitedDataReader::readData( size_t lines_to_skip )
{
    
    std::vector<std::string> tmpChars;
    
    // open file
    std::ifstream readStream( filename.string() );
    if ( not readStream )
    {
        throw RbException()<<"Could not open file "<<filename.make_preferred();
    }
    
    chars.clear();
    
    // read file
    // bool firstLine = true;
    std::string read_line = "";
    size_t lines_skipped = 0;
    while (safeGetline(readStream,read_line))
    {
        ++lines_skipped;
        if ( lines_skipped <= lines_to_skip)
        {
            continue;
        }
        
        // skip blank lines
        std::string::iterator first_nonspace = std::find_if (read_line.begin(), read_line.end(), [](int c) {return not isspace(c);});
        if (first_nonspace == read_line.end())
        {
            continue;
        }

        StringUtilities::stringSplit(read_line, delimiter, tmpChars, true);

        chars.push_back(tmpChars);
        tmpChars.clear();
    };
    
    readStream.close();
}

const std::vector<std::vector<std::string> >& DelimitedDataReader::getChars(void)
{
    return chars;
}


const path& DelimitedDataReader::getFilename(void)
{
    return filename;
}
