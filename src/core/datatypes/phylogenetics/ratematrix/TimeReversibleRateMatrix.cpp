#include <cstddef>
#include <vector>

#include "RbException.h"
#include "TimeReversibleRateMatrix.h"
#include "AbstractRateMatrix.h"
#include "Assignable.h"
#include "Cloneable.h"
#include "MatrixReal.h"
#include "RbVector.h"
#include "RbVectorImpl.h"

using namespace RevBayesCore;

TimeReversibleRateMatrix::TimeReversibleRateMatrix(size_t n) : AbstractRateMatrix(n),
    exchangeability_rates( std::vector<double>(num_states*(num_states-1)/2,1.0) ),
    stationary_freqs( std::vector<double>(num_states,1.0/n) )
{
    
}



TimeReversibleRateMatrix::~TimeReversibleRateMatrix(void)
{
    // nothing to do
}


/**
 * Assign the value of m to this instance. This function is our mechanism to call the assignment operator.
 *
 *
 */
TimeReversibleRateMatrix& TimeReversibleRateMatrix::assign(const Assignable &m)
{
    
    const TimeReversibleRateMatrix *rm = dynamic_cast<const TimeReversibleRateMatrix*>(&m);
    if ( rm != NULL )
    {
        return operator=(*rm);
    }
    else
    {
        throw RbException("Could not assign rate matrix.");
    }
}

/** Calculate the average rate for the rate matrix */
double TimeReversibleRateMatrix::averageRate(void) const
{
    
    double ave = 0.0;
    for (size_t i=0; i<num_states; ++i)
    {
        ave += -stationary_freqs[i] * (*the_rate_matrix)[i][i];
    }
    return ave;
}


void TimeReversibleRateMatrix::computeOffDiagonal( void )
{
    
    MatrixReal& m = *the_rate_matrix;
    
    // set the off-diagonal portions of the rate matrix
    for (size_t i=0, k=0; i<num_states; ++i)
    {
        for (size_t j=i+1; j<num_states; ++j)
        {
            m[i][j] = exchangeability_rates[k] * stationary_freqs[j];
            m[j][i] = exchangeability_rates[k] * stationary_freqs[i];
            k++;
        }
    }
    
    // set flags
    needs_update = true;
}


/** Set the exchangeability rates directly. We assume that we know
 what the exchangeability rates are when this function is called. */
void TimeReversibleRateMatrix::setExchangeabilityRates(const std::vector<double>& er)
{
    
    exchangeability_rates = er;
    
    // set flags
    needs_update = true;
}


/** Set the stationary frequencies directly. We assume that we know
 what the stationary frequencies are when this function is called. */
void TimeReversibleRateMatrix::setStationaryFrequencies(const std::vector<double>& f)
{
    
    stationary_freqs = f;
    
    // set flags
    needs_update = true;
}



std::vector<double> TimeReversibleRateMatrix::getStationaryFrequencies( void ) const
{

    return stationary_freqs;
}

const std::vector<double>& TimeReversibleRateMatrix::getExchangeabilityRates( void ) const
{
    return exchangeability_rates;
}
