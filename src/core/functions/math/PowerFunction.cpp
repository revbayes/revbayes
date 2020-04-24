#include "PowerFunction.h"

#include <cmath>

#include "Cloneable.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;
/**
 *
 *
 * PowerFunction of a TypedDagNode of type double with power of type double.
 *
 *Default Constructor for the Power Function
 * @param b a double for the base of the PowerFunction
 * @param e a double for the exponent of the PowerFcuntion
 */
PowerFunction::PowerFunction(const TypedDagNode<double> *b, const TypedDagNode<double> *e) : ContinuousFunction( new double(0.0) ),
    base( b ),
    exponent( e )
{
    addParameter( b );
    addParameter( e );

}


PowerFunction* PowerFunction::clone( void ) const
{
    return new PowerFunction(*this);
}

void PowerFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == base)
    {
        base = static_cast<const TypedDagNode<double>* >( newP );
    }
    if (oldP == exponent)
    {
        exponent = static_cast<const TypedDagNode<double>* >( newP );
    }
    
}

void PowerFunction::update( void )
{
    *value = pow( base->getValue(), exponent->getValue() );
}


