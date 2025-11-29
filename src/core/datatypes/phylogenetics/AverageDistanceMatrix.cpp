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
DistanceMatrix& AverageDistanceMatrix::getDistanceMatrix(void)
{
    return distanceMatrix;
}


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
MatrixBoolean& AverageDistanceMatrix::getMask(void)
{
    return mask;
}


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


/**
 * Impute missing entries based on the three-point condition satisfied by ultrametric distances
 * (see Hartigan 1967, J. Am. Stat. Assoc. 62(320): 1140--1158, doi:10.1080/01621459.1967.10500922),
 * following De Soete (1984, J. Classif. 1: 235--242, doi:10.1007/BF01890124) and Lapointe & Kirsch
 * (1995, Mol. Biol. Evol. 12(2): 266--284, doi:10.1093/oxfordjournals.molbev.a040209). We perform at
 * most two passes through the matrix, with the distances imputed in the first pass serving as additional
 * input for the second pass.
 */
void AverageDistanceMatrix::ultrametricImputation(void)
{
    int iter = 0;
    int empty = 1;
    
    while (empty > 0)
    {
        // Sanity check: this should never happen
        if (iter > 2)
        {
            throw RbException("More than two iterations are needed to impute all missing entries.");
        }
        
        // Get the row and column indices of all missing entries
        std::vector<size_t> missing_entry_row_indices;
        std::vector<size_t> missing_entry_col_indices;
        
        // Iterate over the lower triangle only
        for (size_t i = 0; i < num_tips; i++)
        {
            for (size_t j = 0; j < i; j++)
            {
                if (mask[i][j] == false)
                {
                    missing_entry_row_indices.push_back(i);
                    missing_entry_col_indices.push_back(j);
                }
            }
        }
        
        /* In each pass, we want to *independently* apply the imputation algorithm to each missing distance
         * d(i, j) for which at least one k can be found such that d(i, k) and d(j, k) are both non-missing.
         * I.e., the imputation of one distance should *not* immediately affect the number of {i, j, k}
         * triplets available for the imputation of the remaining distances. To do so, we create a temporary
         * copy of the matrix, and make element-wise changes to it before re-assigning the original matrix.
         */
        AverageDistanceMatrix* tmp = new AverageDistanceMatrix(*this);
        
        for (size_t i = 0; i < missing_entry_row_indices.size(); i++)
        {
            // For every missing distance d(i, j), find all k such that d(i, k) and d(j, k) are non-missing
            size_t ix0 = missing_entry_row_indices[i];
            size_t ix1 = missing_entry_col_indices[i];
            
            std::vector<size_t> nonmissing;
            
            for (size_t j = 0; j < num_tips; j++)
            {
                if (mask[ix0][j] && mask[ix1][j])
                {
                    nonmissing.push_back(j);
                }
            }
            
            if (nonmissing.size() > 0)
            {
                // For each k, get max[d(i, k), d(j, k)]
                std::vector<double> maxima;
                for (size_t k = 0; k < nonmissing.size(); k++)
                {
                    maxima.push_back( std::max(distanceMatrix[ix0][ nonmissing[k] ], distanceMatrix[ix1][ nonmissing[k] ]) );
                }
                
                // Now, get the minimum over all k
                double repl_dist = *std::min_element(maxima.begin(), maxima.end());
                
                // Re-assign (while exploiting matrix symmetry)
                tmp->getDistanceMatrix()[ix0][ix1] = repl_dist;
                tmp->getDistanceMatrix()[ix1][ix0] = repl_dist;
                tmp->getMask()[ix0][ix1] = true;
                tmp->getMask()[ix1][ix0] = true;
            }
        }
        
        // Re-assign
        distanceMatrix = tmp->getDistanceMatrix();
        mask = tmp->getMask();
        
        // Recalculate the number of missing entries
        empty = (int)missing_entry_row_indices.size();
        
        // Increment counter
        iter++;
    }
    
    std::cout << "Imputation required " << iter - 1 << " iteration(s)." << std::endl;
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
