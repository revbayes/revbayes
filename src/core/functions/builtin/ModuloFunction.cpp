#include "ModuloFunction.h"
#include "Cloneable.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;
/* ModuloFunction of TypedDagNodes of types long
 *
 * @param l a value of type long, the dividend of the remainder division
 * @param r a value of type long, the divisor of the remainder division
 */
ModuloFunction::ModuloFunction(const TypedDagNode<long> *l, const TypedDagNode<long> *r) : TypedFunction<long>( new long(0) ),
    left( l ),
    right( r )
{
    // add the parameters as parents
    this->addParameter( left );
    this->addParameter( right );
    
    update();
}


ModuloFunction::~ModuloFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



ModuloFunction* ModuloFunction::clone( void ) const
{
    return new ModuloFunction( *this );
}


void ModuloFunction::update( void )
{
    
    *this->value = (left->getValue() % right->getValue());
    
}



void ModuloFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if ( oldP == left )
    {
        left = static_cast<const TypedDagNode<long>* >( newP );
    }
    
    if ( oldP == right )
    {
        right = static_cast<const TypedDagNode<long>* >( newP );
    }
    
}

