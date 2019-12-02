#include "DelimitedDataReader.h"

#include <functional>
#include <algorithm>
#include <string>
#include <sstream> // IWYU pragma: keep

#include "RbFileManager.h"
#include "RbException.h"

using namespace RevBayesCore;

DelimitedDataReader::DelimitedDataReader(const std::string &fn, char d, size_t lines_skipped) :
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
    std::ifstream readStream;
    RbFileManager f = RbFileManager(filename);
    if ( f.openFile(readStream) == false )
    {
        throw RbException( "Could not open file " + filename );
    }
    
    chars.clear();
    
    // read file
    // bool firstLine = true;
    std::string read_line = "";
    size_t lines_skipped = 0;
    while (f.safeGetline(readStream,read_line))
    {
        ++lines_skipped;
        if ( lines_skipped <= lines_to_skip)
        {
            continue;
        }
        
        // skip blank lines
        std::string::iterator first_nonspace = std::find_if (read_line.begin(), read_line.end(), std::not1(std::ptr_fun<int,int>(isspace)));
        if (first_nonspace == read_line.end())
        {
            continue;
        }

        std::string field = "";
        std::stringstream ss(read_line);

        int pos = 0;
        while ( std::getline(ss,field,delimiter) )
        {
            field.erase(std::find_if (field.rbegin(), field.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), field.end());
            field.erase(field.begin(), std::find_if (field.begin(), field.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
            tmpChars.push_back(field);
            pos++;
        };
        chars.push_back(tmpChars);
        tmpChars.clear();
    };
    
    f.closeFile( readStream );
}

const std::vector<std::vector<std::string> >& DelimitedDataReader::getChars(void)
{
    return chars;
}


const std::string& DelimitedDataReader::getFilename(void)
{
    return filename;
}
