#include <math.h>
#include <cstddef>
#include <vector>

#include "RateMatrix_F81.h"
#include "TransitionProbabilityMatrix.h"
#include "Assignable.h"
#include "RbException.h"
#include "TimeReversibleRateMatrix.h"

using namespace RevBayesCore;

/** Default Constructor for Felsenstein81 rate matrix
 * @param n the size of the matrix
 *  */
RateMatrix_F81::RateMatrix_F81(size_t n) : TimeReversibleRateMatrix( n )
{
    
    update();
}


/** Destructor */
RateMatrix_F81::~RateMatrix_F81(void)
{
    
}


/**
 * Assign the value of m to this instance. This function is our mechanism to call the assignment operator.
 *
 *
 */
RateMatrix_F81& RateMatrix_F81::assign(const Assignable &m)
{
    
    const RateMatrix_F81 *rm = dynamic_cast<const RateMatrix_F81*>(&m);
    if ( rm != NULL )
    {
        return operator=(*rm);
    }
    else
    {
        throw RbException("Could not assign rate matrix.");
    }
}



/** Calculate the transition probabilities */
void RateMatrix_F81::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const
{
    
    // compute auxilliary variables
    double t = rate * (startAge - endAge);
    double tmp = 1.0;
	for (size_t i=0; i<num_states; i++) tmp -= stationary_freqs[i]*stationary_freqs[i];
    double beta = 1.0 / tmp; 
    double xx = -beta * t;
    double e = exp( xx );
    double oneminuse = 1.0 - e;
    
    // calculate the transition probabilities
	for (size_t i=0; i<num_states; i++) 
    {
        double to_i = stationary_freqs[i] * oneminuse;
		for (size_t j=0; j<num_states; j++) 
        {
            P[j][i] = to_i;    // it is easier to overwrite the case i -> i later than checking for every j
        }
         
        // from i to i
        P[i][i] = stationary_freqs[i] + (1.0-stationary_freqs[i])*e;
    }
}


RateMatrix_F81* RateMatrix_F81::clone( void ) const
{
    return new RateMatrix_F81( *this );
}


void RateMatrix_F81::update( void )
{
    
    if ( needs_update ) 
    {
        // @todo: This is only needed for printing the values of the rate matrix properly to the screen. We should do this more efficiently (Sebastian).
        // We could instead only update the matrix if a print call happened and the matrix was flagged as dirty.
        
        // compute the off-diagonal values
        computeOffDiagonal();
        
        // set the diagonal values
        setDiagonal();
        
        // rescale 
        rescaleToAverageRate( 1.0 );
        
        // clean flags
        needs_update = false;
    }
    
}





