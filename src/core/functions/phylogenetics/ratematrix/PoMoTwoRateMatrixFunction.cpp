#include "PoMoTwoRateMatrixFunction.h"
#include "RateMatrix_PoMoTwo.h"
#include "RbException.h"

using namespace RevBayesCore;


PoMoTwoRateMatrixFunction::PoMoTwoRateMatrixFunction( const TypedDagNode< long > *ps, 
                                                      const TypedDagNode< RbVector<double> > *rho, 
                                                      const TypedDagNode< Simplex > *pi ) : TypedFunction<RateGenerator>( new RateMatrix_PoMoTwo() ), 
populationSize( ps ), 
exchangeabilities( rho ), 
equilibriumFrequencies( pi )

{
   // add the lambda parameter as a parent
   addParameter( populationSize );
   addParameter( exchangeabilities );
   addParameter( equilibriumFrequencies );

   update();
}


PoMoTwoRateMatrixFunction::~PoMoTwoRateMatrixFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


PoMoTwoRateMatrixFunction* PoMoTwoRateMatrixFunction::clone( void ) const
{
    return new PoMoTwoRateMatrixFunction( *this );
}


void PoMoTwoRateMatrixFunction::update( void )
{   
    double p = populationSize ->getValue();
    const std::vector<double>& e = exchangeabilities->getValue();
    const std::vector<double>& f = equilibriumFrequencies->getValue();

    // set parameters
    static_cast< RateMatrix_PoMoTwo* >(value)->setPopulationSize( p );
    static_cast< RateMatrix_PoMoTwo* >(value)->setExchangeabilities( e );
    static_cast< RateMatrix_PoMoTwo* >(value)->setEquilibriumFrequencies( f );
    static_cast< RateMatrix_PoMoTwo* >(value)->update();
}


void PoMoTwoRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == populationSize )
    {
        populationSize = static_cast< const TypedDagNode< long >* >( newP );
    }
    else if (oldP == exchangeabilities)
    {
        exchangeabilities = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == equilibriumFrequencies)
    {
        equilibriumFrequencies = static_cast<const TypedDagNode< Simplex >* >( newP );
    }

}
