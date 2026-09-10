#include "ReversiblePoMoRateMatrixFunction.h"

#include "Cloneable.h"
#include "RateMatrix_ReversiblePoMo.h"
#include "RbVector.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class Simplex; }

using namespace RevBayesCore;


ReversiblePoMoRateMatrixFunction::ReversiblePoMoRateMatrixFunction(const TypedDagNode< std::int64_t > *ps, const TypedDagNode< RbVector<double> > *rho, const TypedDagNode< Simplex > *pi   ) : TypedFunction<RateGenerator>( new RateMatrix_ReversiblePoMo(4 + 6*(ps->getValue() - 1), rho->getValue(), pi->getValue(), ps->getValue() ) ),
    population_size( ps ),
    exchangeabilities( rho ),
    equilibrium_frequencies( pi )
{

    // add the lambda parameter as a parent
    addParameter( population_size );
    addParameter( exchangeabilities );
    addParameter( equilibrium_frequencies );

    update();
}



ReversiblePoMoRateMatrixFunction::~ReversiblePoMoRateMatrixFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



ReversiblePoMoRateMatrixFunction* ReversiblePoMoRateMatrixFunction::clone( void ) const
{
    return new ReversiblePoMoRateMatrixFunction( *this );
}


void ReversiblePoMoRateMatrixFunction::update( void )
{

    // set the exchangeability rates
    static_cast< RateMatrix_ReversiblePoMo* >(value)->setExchangeabilityRates( static_cast< const RbVector<double>& > ( exchangeabilities->getValue() ) );

    // set the stationary frequencies
    static_cast< RateMatrix_ReversiblePoMo* >(value)->setStationaryFrequencies( static_cast< const Simplex& > ( equilibrium_frequencies->getValue() ) );

    // now update the matrix
    static_cast< RateMatrix_ReversiblePoMo* >(value)->update();

}



void ReversiblePoMoRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

    // if (oldP == mutationMatrix)
    // {
    //     mutationMatrix = static_cast<const TypedDagNode< RateGenerator >* >( newP );
    // }
    if (oldP == exchangeabilities)
    {
        exchangeabilities = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == equilibrium_frequencies)
    {
        equilibrium_frequencies = static_cast<const TypedDagNode< Simplex >* >( newP );
    }

}
