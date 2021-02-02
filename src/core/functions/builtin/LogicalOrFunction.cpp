#include "LogicalOrFunction.h"

#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;


/*
 * LogicalOrFunction Constructor for two TypedDagNodes of type Boolean
 *
 * These parameters are the two booleans to be compared by a logical OR
 * @param l A TypedDagNode of type Boolean
 * @param r A TypedDagNode of type Boolean
 */


LogicalOrFunction::LogicalOrFunction(const TypedDagNode<Boolean> *l, const TypedDagNode<Boolean> *r) : TypedFunction<Boolean>( new Boolean(false) ),
    left( l ),
    right( r )
{
    // add the parameters as parents
    this->addParameter( left );
    this->addParameter( right );
    
    update();
}


LogicalOrFunction::~LogicalOrFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



LogicalOrFunction* LogicalOrFunction::clone( void ) const
{
    return new LogicalOrFunction( *this );
}


void LogicalOrFunction::update( void )
{
    
    *this->value = (left->getValue() || right->getValue());
    
}



void LogicalOrFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if ( oldP == left )
    {
        left = static_cast<const TypedDagNode<Boolean>* >( newP );
    }
    if ( oldP == right )
    {
        right = static_cast<const TypedDagNode<Boolean>* >( newP );
    }
    
}

