#include "BinaryMutationCoalescentRateMatrixFunction.h"
#include "RateMatrix_BinaryMutationCoalescent.h"
#include "RbException.h"

using namespace RevBayesCore;

BinaryMutationCoalescentRateMatrixFunction::BinaryMutationCoalescentRateMatrixFunction(size_t n, const TypedDagNode< double > *mr, const TypedDagNode< double > *ne ) : TypedFunction<RateGenerator>( new RateMatrix_BinaryMutationCoalescent(n) ),
    Ne( ne ),
    mu( mr )
{
    
    // add the lambda parameter as a parent
    addParameter( mu );
    addParameter( Ne );
    
    update();
}



BinaryMutationCoalescentRateMatrixFunction::~BinaryMutationCoalescentRateMatrixFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



BinaryMutationCoalescentRateMatrixFunction* BinaryMutationCoalescentRateMatrixFunction::clone( void ) const
{
    return new BinaryMutationCoalescentRateMatrixFunction( *this );
}


void BinaryMutationCoalescentRateMatrixFunction::update( void )
{
    
    // set the base frequencies
    static_cast< RateMatrix_BinaryMutationCoalescent* >(value)->setMutationRate( mu->getValue() );
    static_cast< RateMatrix_BinaryMutationCoalescent* >(value)->setEffectivePopulationSize( Ne->getValue() );
    static_cast< RateMatrix_BinaryMutationCoalescent* >(value)->update();
    
}



void BinaryMutationCoalescentRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == mu)
    {
        mu = static_cast<const TypedDagNode< double >* >( newP );
    }
    else if (oldP == Ne)
    {
        Ne = static_cast<const TypedDagNode< double >* >( newP );
    }
    
}



