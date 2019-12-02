#include "DistanceMatrix.h"

#include <sstream>
#include <string>

#include "DistanceMatrixReader.h"
#include "StringUtilities.h"
#include "RbException.h"
#include "RbVector.h"
#include "RbVectorImpl.h"

using namespace RevBayesCore;

DistanceMatrix::DistanceMatrix( void ) :
    matrix( 2 ),
    taxa( std::vector<Taxon>(2,Taxon()) ),
    num_tips( 2 )
{
    
}

DistanceMatrix::DistanceMatrix( size_t n ) :
    matrix( n ),
    taxa( std::vector<Taxon>(n,Taxon("")) ),
    num_tips( n )
{
    
}


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

const std::vector<Taxon>& DistanceMatrix::getTaxa(void) const
{
    return taxa;
}

const MatrixReal& DistanceMatrix::getMatrix(void) const
{
    return matrix;
}

size_t DistanceMatrix::getSize(void) const
{
	return num_tips;
}


std::string DistanceMatrix::getFilename(void) const
{
    return filename;
}


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


double& DistanceMatrix::getElement( size_t i, size_t j )
{
	return matrix[i][j];
}


void DistanceMatrix::setTaxon(const RevBayesCore::Taxon &t, size_t i)
{
    if ( taxa.size() <= i )
    {
        throw RbException("Cannot set taxon object in distance matrix because of index '"+StringUtilities::to_string(i)+"' out of bounds.");
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
