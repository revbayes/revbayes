#include "LogicalAndFunction.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }
using namespace RevBayesCore;

/*
 * LogicalAndFunction Constructor for two TypedDagNodes of type Boolean
 *
 * These parameters are the two booleans to be compared by a logical AND
 * @param l A TypedDagNode of type Boolean
 * @param r A TypedDagNode of type Boolean
 */

LogicalAndFunction::LogicalAndFunction(const TypedDagNode<Boolean> *l, const TypedDagNode<Boolean> *r) : TypedFunction<Boolean>( new Boolean(false) ),
    left( l ),
    right( r )
{
    // add the parameters as parents
    this->addParameter( left );
    this->addParameter( right );
    
    update();
}


LogicalAndFunction::~LogicalAndFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



LogicalAndFunction* LogicalAndFunction::clone( void ) const
{
    return new LogicalAndFunction( *this );
}


void LogicalAndFunction::update( void )
{
    
    *this->value = (left->getValue() && right->getValue());
    
}



void LogicalAndFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
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

