#include "OrderedRateMatrixFunction.h"
#include "RateMatrix_Ordered.h"
#include "RbException.h"

using namespace RevBayesCore;

OrderedRateMatrixFunction::OrderedRateMatrixFunction(const TypedDagNode<long> *n, const TypedDagNode<double> *l, const TypedDagNode<double> *m) : TypedFunction<RateGenerator>( new RateMatrix_Ordered(n->getValue()) ),
    lambda( l ),
    mu( m )
{
    
    addParameter( lambda );
    addParameter( mu );
    
    update();
}



OrderedRateMatrixFunction::~OrderedRateMatrixFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



OrderedRateMatrixFunction* OrderedRateMatrixFunction::clone( void ) const
{
    return new OrderedRateMatrixFunction( *this );
}


void OrderedRateMatrixFunction::update( void )
{
    double la = lambda->getValue();
    double m = mu->getValue();
    
    static_cast< RateMatrix_Ordered* >(value)->setLambda( la );
    static_cast< RateMatrix_Ordered* >(value)->setMu( m );
    
    static_cast< RateMatrix_Ordered* >(value)->update();
    
}



void OrderedRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == lambda)
    {
        lambda = static_cast<const TypedDagNode<double>* >( newP );
    }
    
    if (oldP == mu)
    {
        mu = static_cast<const TypedDagNode<double>* >( newP );
    }
    
    
}



