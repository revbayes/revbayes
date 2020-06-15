#include "OrderedRateMatrixFunction.h"

#include "RateMatrix_Ordered.h"
#include "Cloneable.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Default constructor.
 *
 * This function takes four inputs:
 * @param n The number of states
 * @param l The rate of gains
 * @param m The number of losses
 * @param allow_zero_state Should state '0' be allowed? (May not be appropriate for some counts)
 */

OrderedRateMatrixFunction::OrderedRateMatrixFunction(const TypedDagNode<long> *n, const TypedDagNode<double> *l, const TypedDagNode<double> *m, bool allow_zero_state) : TypedFunction<RateGenerator>( new RateMatrix_Ordered(n->getValue()) ),
    lambda( l ),
    mu( m )
{
    
    
    
    addParameter( lambda );
    addParameter( mu );
    
    static_cast< RateMatrix_Ordered* >(value)->setAllowZeroState( allow_zero_state );
    
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



