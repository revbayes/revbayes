#include <cstddef>
#include <iosfwd>

#include "AbstractDiscreteTaxonData.h"
#include "PattersonsDFunction.h"
#include "AbstractHomologousDiscreteCharacterData.h"
#include "Cloneable.h"
#include "DiscreteCharacterState.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

PattersonsDFunction::PattersonsDFunction(const std::string &p1, const std::string &p2, const std::string &p3, const std::string &og, const TypedDagNode<AbstractHomologousDiscreteCharacterData> *a) : TypedFunction<double>( new double(0.0) ),
alignment( a ),
P1( p1 ),
P2( p2 ),
P3( p3 ),
OG( og )
{
    // add the lambda parameter as a parent
    addParameter( alignment );
    
    update();
}


PattersonsDFunction::~PattersonsDFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



PattersonsDFunction* PattersonsDFunction::clone( void ) const
{
    
    return new PattersonsDFunction( *this );
}


void PattersonsDFunction::update( void )
{
    size_t num_abba = 0;
    size_t num_baba = 0;
    
    const AbstractDiscreteTaxonData& p1 = alignment->getValue().getTaxonData(P1);
    const AbstractDiscreteTaxonData& p2 = alignment->getValue().getTaxonData(P2);
    const AbstractDiscreteTaxonData& p3 = alignment->getValue().getTaxonData(P3);
    const AbstractDiscreteTaxonData& og = alignment->getValue().getTaxonData(OG);
    
    size_t n = p1.getNumberOfCharacters();
    
    for (size_t i = 0; i < n; i++)
    {
        const DiscreteCharacterState& a = p1[i];
        const DiscreteCharacterState& b = p2[i];
        const DiscreteCharacterState& c = p3[i];
        const DiscreteCharacterState& d = og[i];
        
        // check if this is an abba site
        if (a == d && b == c && a != b)
        {
            ++num_abba;
        }
        
        // check if this is an baba site
        if (a == c && b == d && a != b)
        {
            ++num_baba;
        }
    }
    
    *value = ( double(num_abba) - double(num_baba) ) / double(num_abba + num_baba);
}



void PattersonsDFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == alignment)
    {
        alignment = static_cast<const TypedDagNode< AbstractHomologousDiscreteCharacterData >* >( newP );
    }
    
}
