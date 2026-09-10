#include "MatrixReader.h"

#include <cstdlib>
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


MatrixReader::MatrixReader(const std::string &fn, std::string d, size_t lines_skipped) : DelimitedDataReader(fn, d, lines_skipped)
{
    
    //First, get the size of the matrix
    int nrow = int( chars.size() ) - 1;    //atoi( chars[0][0].c_str() );
    int ncol = int( chars[1].size() ) - 1;
    matrix = MatrixReal( nrow, ncol );

    for (size_t i = 1; i < chars.size(); ++i)
    {
        std::string name = chars[i][0];

        for (size_t j = 1; j < chars[i].size(); ++j)
        {
            matrix[i-1][j-1] = atof( chars[i][j].c_str() );
        }

    }
    
}

const MatrixReal& MatrixReader::getMatrix(void)
{
	return matrix;
}
