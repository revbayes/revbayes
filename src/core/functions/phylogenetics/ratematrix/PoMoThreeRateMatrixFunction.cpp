#include "PoMoThreeRateMatrixFunction.h"
#include "RateMatrix_PoMoThree.h"
#include "RbException.h"

using namespace RevBayesCore;


PoMoThreeRateMatrixFunction::PoMoThreeRateMatrixFunction(const TypedDagNode< long > *ps, const TypedDagNode< RbVector<double> > *rho, const TypedDagNode< Simplex > *pi, const TypedDagNode< RbVector<double> > *gamma   ) : TypedFunction<RateGenerator>( new RateMatrix_PoMoThree() ), 
populationSize( ps ), 
exchangeabilities( rho ), 
equilibriumFrequencies( pi ), 
selectionCoefficients( gamma ) 

{
   // add the lambda parameter as a parent
   addParameter( populationSize );
   addParameter( exchangeabilities );
   addParameter( equilibriumFrequencies );
   addParameter( selectionCoefficients );

   update();
}


PoMoThreeRateMatrixFunction::~PoMoThreeRateMatrixFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


PoMoThreeRateMatrixFunction* PoMoThreeRateMatrixFunction::clone( void ) const
{
    return new PoMoThreeRateMatrixFunction( *this );
}


void PoMoThreeRateMatrixFunction::update( void )
{   
    double p = populationSize ->getValue();
    const std::vector<double>& e = exchangeabilities->getValue();
    const std::vector<double>& f = equilibriumFrequencies->getValue();
    const std::vector<double>& s = selectionCoefficients->getValue();

    // set parameters
    static_cast< RateMatrix_PoMoThree* >(value)->setPopulationSize( p );
    static_cast< RateMatrix_PoMoThree* >(value)->setExchangeabilities( e );
    static_cast< RateMatrix_PoMoThree* >(value)->setEquilibriumFrequencies( f );
    static_cast< RateMatrix_PoMoThree* >(value)->setSelectionCoefficients( s );
    static_cast< RateMatrix_PoMoThree* >(value)->update();
}


void PoMoThreeRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
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
    else if (oldP == selectionCoefficients)
    {
        selectionCoefficients = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
}
