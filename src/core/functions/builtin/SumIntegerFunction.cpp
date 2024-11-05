#include "SumIntegerFunction.h"

#include "RbConstIterator.h"
#include "RbConstIteratorImpl.h"
#include "RbVector.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * SumFunction of a RbVector Constructor.
 * @param v the vector of values of type long
 */
SumIntegerFunction::SumIntegerFunction(const TypedDagNode<RbVector<std::int64_t> > *v) : TypedFunction<std::int64_t>( new long(0.0) ), vals( v )
{
    // add the parameters as parents
    this->addParameter( vals );
    
    update();
}


SumIntegerFunction::~SumIntegerFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



SumIntegerFunction* SumIntegerFunction::clone( void ) const
{
    return new SumIntegerFunction( *this );
}


void SumIntegerFunction::update( void )
{
    
    double m = 0;
    const RbVector<std::int64_t> &v = vals->getValue();
    for ( RbConstIterator<std::int64_t> it = v.begin(); it != v.end(); ++it)
    {
        m += *it;
    }
    
    *this->value = m ;
    
}



void SumIntegerFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if ( oldP == vals )
    {
        vals = static_cast<const TypedDagNode<RbVector<std::int64_t> >* >( newP );
    }
    
}

