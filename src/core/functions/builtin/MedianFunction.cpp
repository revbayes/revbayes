#include "MedianFunction.h"

#include <cstddef>
#include <cmath>

#include "RbVector.h"
#include "RbVectorImpl.h"
#include "TypedDagNode.h"
#include "RbConstants.h" // IWYU pragma: keep

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/** MedianFunction of a RbVector Constructor
 * @param v the vector of values
 */
MedianFunction::MedianFunction(const TypedDagNode< RbVector<double> > *v) : TypedFunction<double>( new double(0.0) ),
    vals( v )
{
    // add the parameters as parents
    this->addParameter( vals );
    
    update();
}


/**
 * Empty destructor.
 */
MedianFunction::~MedianFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



/**
 * Clone function for deep copies.
 */
MedianFunction* MedianFunction::clone( void ) const
{
    return new MedianFunction( *this );
}


/**
 * Update the current value based on the current parameter values.
 */
void MedianFunction::update( void )
{
    
    double m = 0;
    
    // get a copy of the vector which we can sort without changing the original vector
    RbVector<double> v = vals->getValue();

    if ( v.size() == 0 )
    {
        m = RbConstants::Double::nan;
    }
    else if ( v.size() == 1 )
    {
        m = v[0];
    }
    else
    {

        // we need to sort the vector first
        v.sort();

        // do we have an even number of elements
        if (v.size() % 2 == 0)
        {
            size_t index = v.size() / 2;
            m = (v[index-1] + v[index]) / 2.0;
        
        } else {
        
            size_t index = floor(v.size() / 2);
            m = v[index];
        }
        
    }
    
    *this->value = m;
    
}



/**
 * Swap the internal parameters. This happens when the parameters are re-assigned or the entire model graph is cloned.
 * Here we only need to store the new pointer to the vector of real values.
 */
void MedianFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if ( oldP == vals )
    {
        vals = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    
}

