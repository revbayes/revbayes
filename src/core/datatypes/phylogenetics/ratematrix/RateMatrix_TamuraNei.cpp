#include <stddef.h>
#include <cmath>
#include <complex>
#include <vector>

#include "RateMatrix_TamuraNei.h"
#include "RbException.h"
#include "TransitionProbabilityMatrix.h"
#include "Assignable.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "TimeReversibleRateMatrix.h"

using namespace RevBayesCore;

/** Construct rate matrix with 4 states */
RateMatrix_TamuraNei::RateMatrix_TamuraNei(void) : TimeReversibleRateMatrix( 4 )
{
    
    kappa_1 = 1.0;
    kappa_2 = 1.0;
    
    update();
}


/** Destructor */
RateMatrix_TamuraNei::~RateMatrix_TamuraNei(void)
{
    
}


RateMatrix_TamuraNei& RateMatrix_TamuraNei::assign(const Assignable &m)
{
    
    const RateMatrix_TamuraNei *rm = dynamic_cast<const RateMatrix_TamuraNei*>(&m);
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
void RateMatrix_TamuraNei::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const
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
    double beta = 0.5 / (pi_AG * pi_CT + kappa_1 * pi_A * pi_G + kappa_2 * pi_C * pi_T);
    
    // calculate the transition probabilities
    
    double xx = - beta * t;
    double aa = exp(xx);
    double bbR = exp((pi_AG * kappa_1 + pi_CT) * xx);
    double bbY = exp((pi_CT * kappa_2 + pi_AG) * xx);
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


RateMatrix_TamuraNei* RateMatrix_TamuraNei::clone( void ) const
{
    return new RateMatrix_TamuraNei( *this );
}


void RateMatrix_TamuraNei::update( void )
{
    
    if ( needs_update )
    {
        // set the rates (abaaea)
        exchangeability_rates[0] = 1;
        exchangeability_rates[1] = kappa_1;
        exchangeability_rates[2] = 1;
        exchangeability_rates[3] = 1;
        exchangeability_rates[4] = kappa_2;
        exchangeability_rates[5] = 1;
        
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


/**
 * Set the exchangeability rates directly.
 * We assume that we know what the exchangeability rates are when this function is called.
 */
void RateMatrix_TamuraNei::setKappa(double k1, double k2)
{
    
    kappa_1 = k1;
    kappa_2 = k2;
    
    // set flags
    needs_update = true;
}


