#include "SegregatingSitesFunction.h"

#include "AbstractHomologousDiscreteCharacterData.h"
#include "Cloneable.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

SegregatingSitesFunction::SegregatingSitesFunction(const TypedDagNode<AbstractHomologousDiscreteCharacterData> *a, bool excl) : TypedFunction<std::int64_t>( new long(0) ),
    alignment( a ),
    exclude_ambiguous( excl )
{
    // add the lambda parameter as a parent
    addParameter( alignment );
    
    update();
}


SegregatingSitesFunction::~SegregatingSitesFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



SegregatingSitesFunction* SegregatingSitesFunction::clone( void ) const
{
    
    return new SegregatingSitesFunction( *this );
}


void SegregatingSitesFunction::update( void )
{
    
    *value = int( alignment->getValue().getNumberOfSegregatingSites( exclude_ambiguous ) );

}



void SegregatingSitesFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == alignment)
    {
        alignment = static_cast<const TypedDagNode< AbstractHomologousDiscreteCharacterData >* >( newP );
    }
    
}


