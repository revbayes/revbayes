#include "WattersonThetaFunction.h"

#include <cstddef>

#include "RbMathCombinatorialFunctions.h"
#include "AbstractHomologousDiscreteCharacterData.h"
#include "Cloneable.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

WattersonThetaFunction::WattersonThetaFunction(const TypedDagNode<AbstractHomologousDiscreteCharacterData> *a, bool ps, bool excl) : TypedFunction<double>( new double(0.0) ),
    alignment( a ),
    per_site( ps ),
    exclude_ambiguous( excl )
{
    // add the lambda parameter as a parent
    addParameter( alignment );
    
    update();
}


WattersonThetaFunction::~WattersonThetaFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



WattersonThetaFunction* WattersonThetaFunction::clone( void ) const
{
    
    return new WattersonThetaFunction( *this );
}


void WattersonThetaFunction::update( void )
{
    int nss = int( alignment->getValue().getNumberOfSegregatingSites( exclude_ambiguous ) );
    size_t nTaxa = alignment->getValue().getNumberOfTaxa();
    *value = nss / RbMath::harmonicNumber(nTaxa-1);
    
    if ( per_site == true )
    {
        int nSites = int( alignment->getValue().getNumberOfCharacters() );
        *value /= nSites;
    }
    
}



void WattersonThetaFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == alignment)
    {
        alignment = static_cast<const TypedDagNode< AbstractHomologousDiscreteCharacterData >* >( newP );
    }
    
}


