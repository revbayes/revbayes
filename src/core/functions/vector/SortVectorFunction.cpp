#include "SortVectorFunction.h"

#include "TypedDagNode.h"
#include "Cloner.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;


/** Standard constructor from appropriately typed DAG node
 *
 * @param vec a Rbvector of doubles
 * @param a a boolean denoting whether to sort in ascending order
 *
 * */
SortVectorFunction::SortVectorFunction( const TypedDagNode< RbVector<double> >* vec, bool a ) :
    TypedFunction< RbVector<double> >( new RbVector<double>() ),
    ascending( a ),
    realVector( vec )
{
    // add the vector parameter as a parent
    this->addParameter( vec );

    // update the value
    update();
}


/**
 * Destructor. Reference decrement and destruction of parameters is delegated to
 * the deterministic node, so there is nothing to do here.
 */
SortVectorFunction::~SortVectorFunction( void )
{
}


/** Construct type-safe clone of the object. */
SortVectorFunction* SortVectorFunction::clone( void ) const
{
    return new SortVectorFunction( *this );
}


/** Compute the sorted vector from the vector. */
void SortVectorFunction::update( void )
{
    
    // Update unsorted vector
    *value = realVector->getValue();

    value->sort(ascending);
}


/** Swap parameter */
void SortVectorFunction::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{
    if ( oldP == realVector )
        realVector = static_cast< const TypedDagNode< RbVector<double> >* >( newP );
}

