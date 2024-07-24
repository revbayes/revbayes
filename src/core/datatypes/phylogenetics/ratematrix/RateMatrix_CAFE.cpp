#include <math.h>
#include <cstddef>
#include <complex>
#include <vector>

#include "RateMatrix_CAFE.h"
#include "RbMathCombinatorialFunctions.h"
#include "TransitionProbabilityMatrix.h"
#include "MatrixReal.h"
#include "RbException.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "TimeReversibleRateMatrix.h"


using namespace RevBayesCore;

/** Default Constructor for a CAFE rate matrix
 * @param n The size of the matrix
 *
 *  */
RateMatrix_CAFE::RateMatrix_CAFE(size_t n) : TimeReversibleRateMatrix( n+1, false, AbstractRateMatrix::OTHER )
{
    
    binom_coefficients = std::vector< std::vector<double> >( 2*n, std::vector<double>(n, 0) );
    for (size_t i=0; i<(2*n); ++i)
    {
        for (size_t j=0; j<n; ++j)
        {
            binom_coefficients[i][j] = RbMath::choose(i,j);
        }
    }
    
    update();
}



/** Calculate the transition probabilities along a branch
 * @param startAge a double that denotes the start of the branch
 * @param endAge a double that denotes the end of the branch
 * @param rate a double that denotes the overall substitution rate
 * @param P A transition probability matrix
 *
 * Stationary frequencies of nucleotides are also used in this calculation
 *
 *  */
void RateMatrix_CAFE::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const
{
    double t = rate * (startAge - endAge);

    double e = exp((birth-death)*t);
    double alpha = death * (e-1) / (birth * e - death);
    double beta  = birth * (e-1) / (birth * e - death);
    
    bool critical_process = ( birth == death );
    if ( critical_process )
    {
        alpha = birth*t / (1+birth*t);
    }
    for (size_t s=1; s<num_states; ++s)
    {
        double total = 0.0;
        for (size_t c=0; c<num_states; ++c)
        {
            double sum = 0.0;
            
            size_t min = ( s < c ? s : c );
            for (size_t j=0; j <= min; ++j)
            {
                if ( critical_process == true )
                {
                    sum += binom_coefficients[s][j] * binom_coefficients[s+c-j-1][s-1] * pow(alpha, s+c-2*j) * pow(1-2*alpha,j);
                }
                else
                {
                    sum += binom_coefficients[s][j] * binom_coefficients[s+c-j-1][s-1] * pow(alpha, s-j) * pow(beta, c-j) * pow(1-alpha-beta,j);
                }
            }
            
            // sanity check
            if ( sum < 0)
            {
                sum = 0.0;
            }
            
            P[s][c] = sum;
            total += sum;
        }
        
        // now we normalize (because of rounding issues)
        for (size_t c=0; c<num_states; ++c)
        {
            P[s][c] /= total;
        }
        
    }
    
}


RateMatrix_CAFE* RateMatrix_CAFE::clone( void ) const
{
    return new RateMatrix_CAFE( *this );
}



void RateMatrix_CAFE::update( void )
{
    
    if ( needs_update )
    {
        
        // clean flags
        needs_update = false;
    }
    
}


/**
 * Set the rates directly.
 */
void RateMatrix_CAFE::setBirth(double r)
{
    
    birth = r;
    
    // set flags
    needs_update = true;
}


/**
 * Set the rates directly.
 */
void RateMatrix_CAFE::setDeath(double r)
{
    
    death = r;
    
    // set flags
    needs_update = true;
}


