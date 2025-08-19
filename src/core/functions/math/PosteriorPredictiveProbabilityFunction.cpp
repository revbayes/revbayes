#include "PosteriorPredictiveProbabilityFunction.h"

#include <cstddef>

#include "RbVector.h"
#include "RbVectorImpl.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Constructor. Here we simply set up the parameter dependencies.
 */
PosteriorPredictiveProbabilityFunction::PosteriorPredictiveProbabilityFunction(const TypedDagNode< RbVector<double> > *v, const TypedDagNode<double>* k) : TypedFunction<Simplex>( new Simplex(3,1.0) ),
    vals( v ),
    x( k )
{
    // add the parameters as parents
    this->addParameter( vals );
    this->addParameter( x );
    
    update();
}


/**
 * Empty destructor.
 */
PosteriorPredictiveProbabilityFunction::~PosteriorPredictiveProbabilityFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



/**
 * Clone function for deep copies.
 */
PosteriorPredictiveProbabilityFunction* PosteriorPredictiveProbabilityFunction::clone( void ) const
{
    return new PosteriorPredictiveProbabilityFunction( *this );
}


/**
 * Update the current value based on the current parameter values.
 */
void PosteriorPredictiveProbabilityFunction::update( void )
{
    
    const RbVector<double> &v = vals->getValue();
    
    double val = x->getValue();
    size_t size = v.size();
    size_t lower  = 0;
    size_t equal  = 0;
    size_t larger = 0;
    
    for (size_t i=0; i<size; ++i)
    {
        
        if ( v[i] < val )
        {
            ++lower;
        }
        else if ( v[i] == val )
        {
            ++equal;
        }
        else
        {
            ++larger;
        }
        
    }
    
    (*this->value)[0] = lower  / double(size);
    (*this->value)[1] = equal  / double(size);
    (*this->value)[2] = larger / double(size);
    
}



/**
 * Swap the internal parameters. This happens when the parameters are re-assigned or the entire model graph is cloned.
 * Here we only need to store the new pointer to the vector of real values.
 */
void PosteriorPredictiveProbabilityFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if ( oldP == vals )
    {
        vals = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    
    if ( oldP == x )
    {
        x = static_cast<const TypedDagNode<double>* >( newP );
    }
}
