#include "GeneralRateMatrix.h"

#include "RbMathMatrix.h"

using namespace RevBayesCore;


GeneralRateMatrix::GeneralRateMatrix(size_t n, bool rto) : AbstractRateMatrix(n),
//    stationary_freqs( std::vector<double>(num_states,1.0/n) ),
    transition_rates( std::vector<double>(num_states*num_states-num_states, 1.0/n) ),
    rescale_to_one(rto)
{
    
}


GeneralRateMatrix::~GeneralRateMatrix(void)
{
    // nothing to do
}

GeneralRateMatrix& GeneralRateMatrix::assign(const Assignable &m)
{
    const GeneralRateMatrix *rm = dynamic_cast<const GeneralRateMatrix*>(&m);
    if ( rm != NULL )
    {
        return operator=(*rm);
    }
    else
    {
        throw RbException("Could not assign rate matrix.");
    }
}

double GeneralRateMatrix::averageRate(void) const
{
    std::vector<double> stationary_freqs = getStationaryFrequencies();
    
    double ave = 0.0;
    for (size_t i=0; i<num_states; i++)
    {
        ave += -stationary_freqs[i] * (*the_rate_matrix)[i][i];
    }
    
    return ave;
}


std::vector<double> GeneralRateMatrix::getStationaryFrequencies( void ) const
{
    return calculateStationaryFrequencies();
}


/** 
 * Set the exchangeability rates directly. 
 * We assume that we know what the exchangeability rates are when this function is called. 
 */
void GeneralRateMatrix::setTransitionRates(const std::vector<double>& tr)
{
    
    transition_rates = tr;
    
    // set flags
    needs_update = true;
}


void GeneralRateMatrix::update( void )
{
    
    if ( needs_update ) 
    {
        
        // rescale
        if (rescale_to_one)
        {
            rescaleToAverageRate( 1.0 );
        }
        
        // now update the eigensystem
//        updateEigenSystem();
        
        // clean flags
        needs_update = false;
    }
}

