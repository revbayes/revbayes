#include "AverageDistanceMatrix.h"

#include "math.h"
#include <sstream>
#include <string>

#include "StringUtilities.h"
#include "RbException.h"
#include "RbVector.h"
#include "RbVectorImpl.h"

using namespace RevBayesCore;

/** Default constructor for a 2 x 2 average distance matrix */
AverageDistanceMatrix::AverageDistanceMatrix( void ) :
    distanceMatrix( 2 ),
    mask( 2 ),
    num_tips( 2 )
{
    
}

/** Constructor for an n x n average distance matrix */
AverageDistanceMatrix::AverageDistanceMatrix( size_t n ) :
    distanceMatrix( n ),
    mask( n ),
    num_tips( n )
{
    
}


AverageDistanceMatrix::AverageDistanceMatrix(const AverageDistanceMatrix& a)
{
    *this = a;
}

/** Construct an average distance matrix from a distance matrix object and a matrix of Booleans of the same size */
AverageDistanceMatrix::AverageDistanceMatrix(const DistanceMatrix& dm, const MatrixBoolean& m)
{
    distanceMatrix = dm;
    mask = m;
    num_tips = dm.getSize();

}


AverageDistanceMatrix& AverageDistanceMatrix::operator=(const AverageDistanceMatrix& a)
{
    if (this != &a)
    {
        distanceMatrix = a.distanceMatrix;
        mask = a.mask;
        num_tips = a.num_tips;
    }
    
    return *this;
}


AverageDistanceMatrix* AverageDistanceMatrix::clone(void) const
{
    return new AverageDistanceMatrix(*this);
}

/** Get the ratio of defined elements (those marked as 'true' in the Boolean mask) to total elements */
double AverageDistanceMatrix::getCompleteness(void) const
{
    int nonempty = 0;
    for(size_t i = 0; i != num_tips; i++)
    {
        for(size_t j = 0; j != num_tips; j++)
        {
            // Exclude nondefined elements, unless they occur on the diagonal
            if(i == j || mask[i][j]) nonempty++;
        }
    }
    
    double completeness = static_cast<double>(nonempty) / pow(num_tips, 2.0);
    return completeness;
}

/** Get the distance matrix (without the mask distinguishing between defined and undefined entries) */
const DistanceMatrix& AverageDistanceMatrix::getDistanceMatrix(void) const
{
    return distanceMatrix;
}

/** Extract a distance from a matrix and check whether it is defined or not
 * @param i Row index
 * @param j Column index
 *
 * @return A pair of a raw real-valued distance and a Boolean indicating whether it is defined or not
 */
std::pair<double, bool> AverageDistanceMatrix::getElement( size_t i, size_t j )
{
    return std::make_pair( distanceMatrix.getElement(i,j), mask[i][j] );
}

/** Get the Boolean mask distinguishing between defined (true) and undefined (false) elements */
const MatrixBoolean& AverageDistanceMatrix::getMask(void) const
{
    return mask;
}

/** Get the number of tips of the tree associated with the matrix */
size_t AverageDistanceMatrix::getSize(void) const
{
    return num_tips;
}

/** Get the taxa whose pairwise distances are stored in the matrix */
const std::vector<Taxon>& AverageDistanceMatrix::getTaxa(void) const
{
    return distanceMatrix.getTaxa();
}

/** Get the number of elements in a row or column of the matrix */
size_t AverageDistanceMatrix::size(void) const
{
    return distanceMatrix.size();
}


std::ostream& RevBayesCore::operator<<(std::ostream& o, const AverageDistanceMatrix& x)
{
    
    std::stringstream s;
    
    // Generate nice header
    o << std::endl;

    s << "AverageDistanceMatrix with " << x.getSize() << " tips. " << std::endl;

    o << s.str();
    std::vector<Taxon> taxa = x.getTaxa();

    for ( size_t i = 0; i < x.getSize(); ++i )
    {
        o << taxa[i] ;
        for ( size_t j = 0; j < x.getSize(); ++j )
        {
            if (x.getMask()[i][j])
            {
                o << "\t" << x.getDistanceMatrix()[i][j] << " T" ;
            } else {
                o << "\t" << x.getDistanceMatrix()[i][j] << " F" ;
            }
        }
        o << std::endl;
    }
    o << std::endl;
    
    return o;
}
