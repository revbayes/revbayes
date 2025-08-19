#include "DistanceMatrixReader.h"

#include <cstdlib>
#include <string>

#include "RbVector.h"
#include "RbVectorImpl.h"


using namespace RevBayesCore;


DistanceMatrixReader::DistanceMatrixReader(const std::string &fn, std::string d, size_t numSkipped) : DelimitedDataReader(fn, d, numSkipped)
{
	
	//First, get the size of the matrix
	int siz = int(chars.size()) -1;//atoi( chars[0][0].c_str() );
	matrix = MatrixReal( siz );

	for (size_t i = 1; i < chars.size(); ++i)
	{
		std::string name = chars[i][0];
		
		taxa.push_back( Taxon(name) );
		for (size_t j = 1; j < chars[i].size(); ++j)
		{
			matrix[i-1][j-1] = atof( chars[i][j].c_str() );
		}
		
	}
}


const std::vector<Taxon>& DistanceMatrixReader::getTaxa(void)
{
	return taxa;
}


const MatrixReal& DistanceMatrixReader::getMatrix(void)
{
	return matrix;
}
