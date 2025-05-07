#include "EmpiricalQuantileFunction.h"

#include <cstddef>
#include <cmath>

#include "RbVector.h"
#include "RbVectorImpl.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;



/**
 * EmpiricalQuantileFunction of a RbVector constructor and a TypedDagNode
 * double of the kth quantile
 *
 * @param v an RbVector of values of type double
 * @param k a double value
 */
EmpiricalQuantileFunction::EmpiricalQuantileFunction(const TypedDagNode< RbVector<double> > *v, const TypedDagNode<double>* k) : TypedFunction<double>( new double(0.0) ),
    vals( v ),
    kth_quantile( k )
{
    // add the parameters as parents
    this->addParameter( vals );
    this->addParameter( kth_quantile );
    
    update();
}


/**
 * Empty destructor.
 */
EmpiricalQuantileFunction::~EmpiricalQuantileFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



/**
 * Clone function for deep copies.
 */
EmpiricalQuantileFunction* EmpiricalQuantileFunction::clone( void ) const
{
    return new EmpiricalQuantileFunction( *this );
}


/**
 * Update the current value based on the current parameter values.
 */
void EmpiricalQuantileFunction::update( void )
{
    
    RbVector<double> v = vals->getValue();
    v.sort();
//    std::sort(v.begin(),v.end());
    
    std::size_t index = round( (v.size()-1) * kth_quantile->getValue());
    
    *this->value = v[index];
    
}



/**
 * Swap the internal parameters. This happens when the parameters are re-assigned or the entire model graph is cloned.
 * Here we only need to store the new pointer to the vector of real values.
 */
void EmpiricalQuantileFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if ( oldP == vals )
    {
        vals = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    
    if ( oldP == kth_quantile )
    {
        kth_quantile = static_cast<const TypedDagNode<double>* >( newP );
    }
}
