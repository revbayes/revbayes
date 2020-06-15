#include "LogicalUnotFunction.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }
using namespace RevBayesCore;

LogicalUnotFunction::LogicalUnotFunction(const TypedDagNode<Boolean> *l) : TypedFunction<Boolean>( new Boolean(false) ),
    left( l )
{
    // add the parameters as parents
    this->addParameter( left );
    
    update();
}


LogicalUnotFunction::~LogicalUnotFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



LogicalUnotFunction* LogicalUnotFunction::clone( void ) const
{

    return new LogicalUnotFunction( *this );
}


void LogicalUnotFunction::update( void )
{
    
    *(this->value) = !( left->getValue() );
    
}



void LogicalUnotFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {
    
    if ( oldP == left ) {
        left = static_cast<const TypedDagNode<Boolean>* >( newP );
    }
    
}

