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


MatrixReader::MatrixReader(const std::string &fn, std::string d, size_t lines_skipped, bool hasHeaders) : DelimitedDataReader(fn, d, lines_skipped)
{
    
    const size_t row_start = hasHeaders ? 1 : 0;
    const size_t col_start = hasHeaders ? 1 : 0;

    int nrow = int( chars.size() ) - int( row_start );
    int ncol = int( chars[row_start].size() ) - int( col_start );
    matrix = MatrixReal( nrow, ncol );

    for (size_t i = row_start; i < chars.size(); ++i)
    {
        for (size_t j = col_start; j < chars[i].size(); ++j)
        {
            matrix[i - row_start][j - col_start] = atof( chars[i][j].c_str() );
        }
    }
    
}

const MatrixReal& MatrixReader::getMatrix(void)
{
	return matrix;
}
