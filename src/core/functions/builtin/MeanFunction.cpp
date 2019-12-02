#include "MeanFunction.h"

using namespace RevBayesCore;

/** MeanFunction of a RbVector Constructor
 * @param v the vector of values
 */
MeanFunction::MeanFunction(const TypedDagNode< RbVector<double> > *v) : TypedFunction<double>( new double(0.0) ),
    vals( v )
{
    // add the parameters as parents
    this->addParameter( vals );
    
    update();
}


/**
 * Empty destructor.
 */
MeanFunction::~MeanFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



/**
 * Clone function for deep copies.
 */
MeanFunction* MeanFunction::clone( void ) const
{
    return new MeanFunction( *this );
}


/**
 * Update the current value based on the current parameter values.
 */
void MeanFunction::update( void )
{
    
    double m = 0;
    const std::vector<double> &v = vals->getValue();
    for ( std::vector<double>::const_iterator it = v.begin(); it != v.end(); ++it)
    {
        m += *it;
    }
    
    *this->value = m / v.size();
    
}



/**
 * Swap the internal parameters. This happens when the parameters are re-assigned or the entire model graph is cloned.
 * Here we only need to store the new pointer to the vector of real values.
 */
void MeanFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if ( oldP == vals )
    {
        vals = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    
}

