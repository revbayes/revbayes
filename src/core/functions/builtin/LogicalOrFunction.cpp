#include "LogicalOrFunction.h"

using namespace RevBayesCore;

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

