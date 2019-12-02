#include "RateMatrix_JC.h"
#include "RbException.h"
#include "RbMathMatrix.h"
#include "TransitionProbabilityMatrix.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>

using namespace RevBayesCore;

/** Construct rate matrix with n states */
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


/**
 * Assign the value of m to this instance. This function is our mechanism to call the assignment operator.
 *
 *
 */
RateMatrix_JC& RateMatrix_JC::assign(const Assignable &m)
{
    
    const RateMatrix_JC *rm = dynamic_cast<const RateMatrix_JC*>(&m);
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

