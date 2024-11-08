#include "RateMatrix_JC.h"

#include <math.h>

#include "TransitionProbabilityMatrix.h"
#include "RbException.h"





using namespace RevBayesCore;

/** Jukes Cantor Rate Matrix Constructor
 *
 * @param n the size of the matrix
 *  */
RateMatrix_JC::RateMatrix_JC(size_t n) : TimeReversibleRateMatrix( n )
{
    
    
    // compute the off-diagonal values
    computeOffDiagonal();
            
    // set the diagonal values
    setDiagonal();
            
    // rescale 
    rescaleToAverageRate( 1.0 );
            
}


/** Destructor */
RateMatrix_JC::~RateMatrix_JC(void)
{
    
}

/** Calculate the transition probabilities */
void RateMatrix_JC::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const
{
    
    
    double t = rate * (startAge - endAge);
    
    // calculate the transition probabilities
    double bf = 1.0 / num_states;
    double oneMinusBf = 1.0 - bf;
    double p_ii = bf + oneMinusBf * exp(-t/oneMinusBf);
    double p_ij = bf - bf * exp(-t/oneMinusBf);
	for (size_t i=0; i<num_states; i++) 
    {
        P[i][i] = p_ii;
		for (size_t j=i+1; j<num_states; j++) 
        {
            P[i][j] = p_ij;
            P[j][i] = p_ij;
        }
        
    }
    
}


RateMatrix_JC* RateMatrix_JC::clone( void ) const
{
    return new RateMatrix_JC( *this );
}


void RateMatrix_JC::update( void )
{
    // nothing to do
}

