#include "MatrixReader.h"

#include <stdlib.h>
#include <algorithm>
#include <functional>
#include <istream>
#include <string>
#include <vector>

#include "RbFileManager.h"
#include "RbException.h"
#include "RbVector.h"
#include "RbVectorImpl.h"


using namespace RevBayesCore;


MatrixReader::MatrixReader(const std::string &fn, char d, size_t lines_skipped) :
    filename(fn),
    delimiter(d)
{
    
    readData( lines_skipped );
    
}

const std::string& MatrixReader::getFilename(void)
{
    return filename;
}

const MatrixReal& MatrixReader::getMatrix(void)
{
	return matrix;
}

void MatrixReader::readData( size_t lines_to_skip )
{
    
    
    // make a growable container
    std::vector<std::string> tmpChars;
    std::vector<std::vector<double> > tmpMatrix;
    
    // open file
    std::ifstream readStream;
    RbFileManager f = RbFileManager(filename);
    if ( f.openFile(readStream) == false )
    {
        throw RbException( "Could not open file " + filename );
    }
    
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
        }
        
        double nchar = tmpChars.size();
        std::vector<double> tmpValues(nchar);
        for (size_t i = 0; i < nchar; ++i) {
            tmpValues[i] = atof( tmpChars[i].c_str() );
        }
        tmpMatrix.push_back(tmpValues);
        
        tmpChars.clear();
 
    }


    // make sure that all the rows are of the same size
    size_t nrow = tmpMatrix.size();
    size_t ncol = tmpMatrix[0].size();
    for (size_t i = 1; i < nrow; ++i)
    {
        if ( tmpMatrix[i].size() != ncol )
        {
            RbException("Error reading matrix file: not all rows contain the same number of columns.");
        }
    }
    
    
    
    // make the matrix
    matrix = MatrixReal(nrow, ncol);
    for (size_t i = 0; i < nrow; ++i)
    {
        for (size_t j = 0; j < ncol; ++j)
        {
            matrix[i][j] = tmpMatrix[i][j];
        }
    }
    

    // close the input file connection
    f.closeFile( readStream );
    
}
