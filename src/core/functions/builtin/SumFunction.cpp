#include "SumFunction.h"

using namespace RevBayesCore;
/**
 * SumFunction of a RbVector Constructor.
 * @param v the vector of values of type double
 */
SumFunction::SumFunction(const TypedDagNode<RbVector<double> > *v) : TypedFunction<double>( new double(0.0) ), vals( v )
{
    // add the parameters as parents
    this->addParameter( vals );
    
    update();
}


SumFunction::~SumFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



SumFunction* SumFunction::clone( void ) const
{
    return new SumFunction( *this );
}


void SumFunction::update( void )
{
    
    double m = 0;
    const RbVector<double> &v = vals->getValue();
    for ( RbConstIterator<double> it = v.begin(); it != v.end(); ++it)
    {
        m += *it;
    }
    
    *this->value = m ;
    
}



void SumFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if ( oldP == vals )
    {
        vals = static_cast<const TypedDagNode<RbVector<double> >* >( newP );
    }
    
}

