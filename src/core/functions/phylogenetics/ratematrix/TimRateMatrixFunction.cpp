#include "TimRateMatrixFunction.h"
#include "RbException.h"

using namespace RevBayesCore;

TimRateMatrixFunction::TimRateMatrixFunction(const TypedDagNode< Simplex > *er, const TypedDagNode< Simplex > *bf) : TypedFunction<RateGenerator>( new RateMatrix_TIM(bf->getValue().size()) ),
    exchangeability_rates( er ),
    base_frequencies( bf )
{
    // add the lambda parameter as a parent
    addParameter( base_frequencies );
    addParameter( exchangeability_rates );
    
    update();
}

TimRateMatrixFunction::~TimRateMatrixFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



TimRateMatrixFunction* TimRateMatrixFunction::clone( void ) const
{
    return new TimRateMatrixFunction( *this );
}


void TimRateMatrixFunction::update( void )
{
    
    // get the information from the arguments for reading the file
    const std::vector<double>& r = exchangeability_rates->getValue();
    const std::vector<double>& f = base_frequencies->getValue();
    
    
    // set the base frequencies
    static_cast< RateMatrix_TIM* >(value)->setStationaryFrequencies( f );
    static_cast< RateMatrix_TIM* >(value)->setRates( r );
    
    value->update();
    
}



void TimRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == base_frequencies)
    {
        base_frequencies = static_cast<const TypedDagNode< Simplex >* >( newP );
    }
    else if (oldP == exchangeability_rates)
    {
        exchangeability_rates = static_cast<const TypedDagNode< Simplex >* >( newP );
    }
    
}


