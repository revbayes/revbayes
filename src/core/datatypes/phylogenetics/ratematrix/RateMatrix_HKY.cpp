#include <cstddef>
#include <cmath>
#include <vector>

#include "MatrixReal.h"
#include "RateMatrix_HKY.h"
#include "RbException.h"
#include "TransitionProbabilityMatrix.h"
#include "Assignable.h"
#include "Cloneable.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "TimeReversibleRateMatrix.h"

using namespace RevBayesCore;

/** Default Constructor for a HKY rate matrix
 * Kappa is initialized to 1
 *  */
RateMatrix_HKY::RateMatrix_HKY(void) : TimeReversibleRateMatrix( 4 )
{
    
    kappa = 1.0;
    
    update();
    
}


/** Destructor */
RateMatrix_HKY::~RateMatrix_HKY(void)
{
    
}


/**
 * Assign the value of m to this instance. This function is our mechanism to call the assignment operator.
 *
 *
 */
RateMatrix_HKY& RateMatrix_HKY::assign(const Assignable &m)
{
    
    const RateMatrix_HKY *rm = dynamic_cast<const RateMatrix_HKY*>(&m);
    if ( rm != NULL )
    {
        return operator=(*rm);
    }
    else
    {
        throw RbException("Could not assign rate matrix.");
    }
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
void RateMatrix_HKY::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const
{
    
    double t = rate * (startAge - endAge);
    
    // notation:
    double pi_A = stationary_freqs[0];
    double pi_C = stationary_freqs[1];
    double pi_G = stationary_freqs[2];
    double pi_T = stationary_freqs[3];
    
    // compute auxilliary variables
    double pi_AG = pi_A + pi_G;
    double pi_CT = pi_C + pi_T;
    
    // compute beta
    double beta = 1.0 / (2.0*pi_AG*pi_CT+2.0*kappa*((pi_A*pi_G)+(pi_C*pi_T))); 
    
    // calculate the transition probabilities
    
    double xx = - beta * t;
    double aa = exp(xx);
    double bbR = exp( (1.0+pi_AG*(kappa-1.0))* xx);
    double bbY = exp( (1.0+pi_CT*(kappa-1.0))* xx);
    double oneminusa = 1 - aa;
    
    P[0][0] = (pi_A * ( pi_AG+pi_CT*aa ) + pi_G*bbR) / pi_AG;
    P[0][1] = pi_C * oneminusa;
    P[0][2] = (pi_G * ( pi_AG+pi_CT*aa ) - pi_G*bbR) / pi_AG;
    P[0][3] = pi_T * oneminusa;
    
    P[1][0] = pi_A * oneminusa;
    P[1][1] = (pi_C * ( pi_CT+pi_AG*aa ) + pi_T*bbY) / pi_CT;
    P[1][2] = pi_G * oneminusa;
    P[1][3] = (pi_T * ( pi_CT+pi_AG*aa ) - pi_T*bbY) / pi_CT;
    
    P[2][0] = (pi_A * ( pi_AG+pi_CT*aa ) - pi_A*bbR) / pi_AG;
    P[2][1] = P[0][1];
    P[2][2] = (pi_G * ( pi_AG+pi_CT*aa ) + pi_A*bbR) / pi_AG;
    P[2][3] = P[0][3];
    
    P[3][0] = P[1][0];
    P[3][1] = (pi_C * ( pi_CT+pi_AG*aa ) - pi_C*bbY) / pi_CT;
    P[3][2] = P[1][2];
    P[3][3] = (pi_T * ( pi_CT+pi_AG*aa ) + pi_C*bbY) / pi_CT;
	
}


RateMatrix_HKY* RateMatrix_HKY::clone( void ) const
{
    
    return new RateMatrix_HKY( *this );
}


void RateMatrix_HKY::setKappa( double k )
{
    
    kappa = k;
    
    // set flags
    needs_update = true;
    
}


void RateMatrix_HKY::update( void )
{
    
    if ( needs_update ) 
    {
        MatrixReal &m = *the_rate_matrix;
        
        // @todo: This is only needed for printing the values of the rate matrix properly to the screen. We should do this more efficiently (Sebastian).
        // We could instead only update the matrix if a print call happened and the matrix was flagged as dirty.
        
        // compute the off-diagonal values
        m[0][1] = stationary_freqs[1];
        m[0][2] = kappa*stationary_freqs[2];
        m[0][3] = stationary_freqs[3];
        
        m[1][0] = stationary_freqs[0];
        m[1][2] = stationary_freqs[2];
        m[1][3] = kappa*stationary_freqs[3];
        
        m[2][0] = kappa*stationary_freqs[0];
        m[2][1] = stationary_freqs[1];
        m[2][3] = stationary_freqs[3];
        
        m[3][0] = stationary_freqs[0];
        m[3][1] = kappa*stationary_freqs[1];
        m[3][2] = stationary_freqs[2];
        
        // set the diagonal values
        setDiagonal();
        
        // rescale 
        rescaleToAverageRate( 1.0 );
        
        // clean flags
        needs_update = false;
    }
    
}


