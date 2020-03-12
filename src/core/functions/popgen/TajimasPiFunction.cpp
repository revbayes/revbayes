#include "TajimasPiFunction.h"

#include "AbstractHomologousDiscreteCharacterData.h"
#include "Cloneable.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/*
 * Tajimas Pi Function Constructor
 *
 * @param a A character data set with the alignment
 * @param ps A boolean for whether we calculate Tajimas Pi per site
 * @param excl A boolean for whether we exclude the ambiguous sites from the alignment
 */


TajimasPiFunction::TajimasPiFunction(const TypedDagNode<AbstractHomologousDiscreteCharacterData> *a, bool ps, bool excl) : TypedFunction<double>( new double(0.0) ),
    alignment( a ),
    per_site( ps ),
    exclude_ambiguous_sites( excl )
{
    // add the lambda parameter as a parent
    addParameter( alignment );
    
    update();
}


TajimasPiFunction::~TajimasPiFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



TajimasPiFunction* TajimasPiFunction::clone( void ) const
{
    
    return new TajimasPiFunction( *this );
}


void TajimasPiFunction::update( void )
{
    double pd  = alignment->getValue().getAveragePaiwiseSequenceDifference( exclude_ambiguous_sites );
    *value = pd;
    
    if ( per_site == true )
    {
        int nSites = int( alignment->getValue().getNumberOfCharacters() );
        *value /= nSites;
    }
    
}



void TajimasPiFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == alignment)
    {
        alignment = static_cast<const TypedDagNode< AbstractHomologousDiscreteCharacterData >* >( newP );
    }
    
}


