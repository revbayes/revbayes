#include "RateMatrix_InfiniteSites.h"

#include <cmath>

#include "RbException.h"
#include "TransitionProbabilityMatrix.h"
#include "Cloneable.h"
#include "MatrixReal.h"
#include "RbVector.h"
#include "RbVectorImpl.h"

using namespace RevBayesCore;

/** Construct rate matrix with n states */
RateMatrix_InfiniteSites::RateMatrix_InfiniteSites(size_t n) : TimeReversibleRateMatrix( n )
{
    
    MatrixReal &m = *the_rate_matrix;
    
    // compute the off-diagonal values
    for (size_t i=1; i<num_states; i++)
    {
        m[0][i] = 1.0;
        for (size_t j=0; j<num_states; j++)
        {
            m[i][j] = 0;
        }
    }
    
    // set the diagonal values
    setDiagonal();
    
    // rescale
    rescaleToAverageRate( 1.0 );
    
}


/** Destructor */
RateMatrix_InfiniteSites::~RateMatrix_InfiniteSites(void)
{
    
}


/** Calculate the transition probabilities */
void RateMatrix_InfiniteSites::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const
{
    
    double t = rate * (startAge - endAge);
    
    // calculate the transition probabilities
    P[0][0] = exp(-t);
    for (size_t i=1; i<num_states; i++)
    {
        P[0][i] = (1 - exp(-t)) / (num_states - 1.0);
        for (size_t j=0; j<num_states; j++)
        {
            P[i][j] = 0;
        }
        P[i][i] = 1;
    }
}


RateMatrix_InfiniteSites* RateMatrix_InfiniteSites::clone( void ) const
{
    return new RateMatrix_InfiniteSites( *this );
}


void RateMatrix_InfiniteSites::update( void )
{
    // nothing to do
}

