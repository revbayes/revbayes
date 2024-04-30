#include "GeneralRateMatrix.h"

#include "Cloneable.h"
#include "MatrixReal.h"
#include "RbException.h"
#include "RbVector.h"
#include "RbVectorImpl.h"

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

double GeneralRateMatrix::averageRate(void) const
{
    std::vector<double> stationary_freqs = getStationaryFrequencies();

    std::vector<int> emit = get_emitted_letters();
    
    double ave = 0.0;
    for (size_t i=0; i<num_states; i++)
    {
        for (size_t j=0; j<num_states; j++)
        {
            if (emit[i] != emit[j])
            {
                ave += stationary_freqs[i] * (*the_rate_matrix)[i][j];
            }
            
        }
        
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

