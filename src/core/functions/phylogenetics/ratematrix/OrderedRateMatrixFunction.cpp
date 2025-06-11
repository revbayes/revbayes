#include "OrderedRateMatrixFunction.h"

#include "RateMatrix_FreeK.h"
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

OrderedRateMatrixFunction::OrderedRateMatrixFunction(const TypedDagNode<std::int64_t> *n, const TypedDagNode<double> *l, const TypedDagNode<double> *m, bool allow_zero_state, bool rescale, std::string method) : TypedFunction<RateGenerator>( new RateMatrix_FreeK( n->getValue(), rescale, method ) ),
    lambda( l ),
    mu( m ),
    zero( allow_zero_state )
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
    double l = lambda->getValue();
    double m = mu->getValue();

    size_t n = static_cast< RateMatrix_FreeK* >(value)->getNumberOfStates();
    
    std::vector<double> r_flat( n * (n-1) );
    size_t k = 0;

    for (size_t i=0; i< n; i++)
    {
        for (size_t j=0; j< n; j++)
        {
            if ( zero == true || (j != 0 && i != 0) )
            {
                if (j == i+1)
                {
                    r_flat[k] = l;
                }
                else if (j == i-1)
                {
                    r_flat[k] = m;
                }
            }

            if (j != i)
            {
                k++;
            }
        }
    }

    // set the flattened rates
    static_cast< RateMatrix_FreeK* >(value)->setTransitionRates(r_flat);

    value->update();
    
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



