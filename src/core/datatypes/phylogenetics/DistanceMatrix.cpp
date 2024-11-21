#include "DistanceMatrix.h"

#include <sstream>
#include <string>

#include "DistanceMatrixReader.h"
#include "StringUtilities.h"
#include "RbException.h"
#include "RbVector.h"
#include "RbVectorImpl.h"

using namespace RevBayesCore;

/** Default constructor for a 2 x 2 distance matrix */
DistanceMatrix::DistanceMatrix( void ) :
    matrix( 2 ),
    taxa( std::vector<Taxon>(2,Taxon()) ),
    num_tips( 2 )
{
    
}

/** Constructor for an n x n distance matrix */
DistanceMatrix::DistanceMatrix( size_t n ) :
    matrix( n ),
    taxa( std::vector<Taxon>(n,Taxon("")) ),
    num_tips( n )
{
    
}

/** Read in a distance matrix from a file */
DistanceMatrix::DistanceMatrix(DistanceMatrixReader* tadr) : filename(tadr->getFilename())
{
    taxa = tadr->getTaxa();
	matrix = tadr->getMatrix();
	num_tips = taxa.size();
}

DistanceMatrix::DistanceMatrix(const DistanceMatrix& a)
{
    *this = a;
}

/** Construct a distance matrix from a real-valued matrix and a vector of taxa of the corresponding size */
DistanceMatrix::DistanceMatrix(const MatrixReal& a, const std::vector<Taxon>& t)
{
	taxa = t;
	matrix = a;
	num_tips = taxa.size();

}

DistanceMatrix& DistanceMatrix::operator=(const DistanceMatrix& a)
{
    if (this != &a)
    {
        taxa = a.taxa;
		matrix = a.matrix;
		filename = a.filename;
		num_tips = a.num_tips;
    }
    
    return *this;
}

DistanceMatrix* DistanceMatrix::clone(void) const
{
    return new DistanceMatrix(*this);
}

/** Get the taxa whose pairwise distances are stored in the matrix */
const std::vector<Taxon>& DistanceMatrix::getTaxa(void) const
{
    return taxa;
}

/** Get the real-valued matrix of distances */
const MatrixReal& DistanceMatrix::getMatrix(void) const
{
    return matrix;
}

/** Get the number of tips of the tree associated with the matrix */
size_t DistanceMatrix::getSize(void) const
{
	return num_tips;
}


const path& DistanceMatrix::getFilename(void) const
{
    return filename;
}

/** Get the number of elements in a row or column of the matrix */
size_t DistanceMatrix::size(void) const
{
	return matrix.size();
}


RbVector<double>& DistanceMatrix::operator[]( size_t index )
{
	
	return matrix[index];
}


const RbVector<double>& DistanceMatrix::operator[]( size_t index ) const
{
	return matrix[index];
}

/** Extract a distance from a matrix
 * @param i Row index
 * @param j Column index
 */
double& DistanceMatrix::getElement( size_t i, size_t j )
{
	return matrix[i][j];
}

/** Set a taxon
 * @param t Taxon that will be set if index i is not out of bound
 * @param i Index denoting the element of the 'taxa' vector that is to be set
 *
 * @throw RbException if i is out of bounds
 */
void DistanceMatrix::setTaxon(const RevBayesCore::Taxon &t, size_t i)
{
    if ( taxa.size() <= i )
    {
        throw RbException() << "Cannot set taxon object in distance matrix because of index '" << i << "' out of bounds."; 
    }
    taxa[i] = t;
}


std::ostream& RevBayesCore::operator<<(std::ostream& o, const DistanceMatrix& x)
{
	
    std::stringstream s;
    
    // Generate nice header
    o << std::endl;

    s << "DistanceMatrix with " << x.getSize() << " tips. " << std::endl;

    o << s.str();
	std::vector<Taxon> taxa = x.getTaxa();

	for ( size_t i = 0; i < x.getSize(); ++i )
    {
		o << taxa[i] ;
		for ( size_t j = 0; j < x.getSize(); ++j )
        {
        	o << "\t" << x.getMatrix()[i][j] ;
		}
		o << std::endl;
	}
    o << std::endl;
    
    return o;
}
