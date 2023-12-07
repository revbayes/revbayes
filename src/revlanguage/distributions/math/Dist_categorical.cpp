#include <math.h>
#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "CategoricalDistribution.h"
#include "Dist_categorical.h"
#include "RlSimplex.h"
#include "DagMemberFunction.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "Natural.h"
#include "RbHelpReference.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlDagMemberFunction.h"
#include "RlDeterministicNode.h"
#include "RlDistribution.h"
#include "RlStochasticNode.h"
#include "RlTypedDistribution.h"
#include "RlTypedFunction.h"
#include "StochasticNode.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "TypedFunction.h"

namespace RevBayesCore { class Simplex; }

using namespace RevLanguage;

Dist_categorical::Dist_categorical() :
    TypedDistribution<Natural>()
{
    setGuiDistributionName("Categorical");
    setGuiDistributionToolTip("Categorical distribution for random variables on the space {1, 2, ..., k}");
}


Dist_categorical::~Dist_categorical()
{
    
}



Dist_categorical* Dist_categorical::clone( void ) const
{
    return new Dist_categorical(*this);
}


RevBayesCore::CategoricalDistribution* Dist_categorical::createDistribution( void ) const
{
    // get the parameters
    RevBayesCore::TypedDagNode< RevBayesCore::Simplex >*            q    = static_cast<const Simplex&>( p->getRevObject() ).getDagNode();
    RevBayesCore::CategoricalDistribution*                          d    = new RevBayesCore::CategoricalDistribution( q );
    
    return d;
}


/* Get Rev type of class (static) */
const std::string& Dist_categorical::getClassType( void )
{
    static std::string rev_type = "Dist_categorical";
    
    return rev_type;
}


/* Get Rev type spec of class (static) */
const TypeSpec& Dist_categorical::getClassTypeSpec(void)
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Distribution::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_categorical::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "cat" );
    
    return a_names;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_categorical::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "categorical";
    
    return d_name;
}


/** Return member rules */
const MemberRules& Dist_categorical::getParameterRules(void) const
{
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        memberRules.push_back( new ArgumentRule( "p", Simplex::getClassTypeSpec(), "The probability for each category.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return memberRules;
}


/** Get Rev type spec of instance (dynamic) */
const TypeSpec& Dist_categorical::getTypeSpec( void ) const
{
    return getClassTypeSpec();
}


/** Print value for user */
void Dist_categorical::printValue(std::ostream& o) const
{
    
    o << "Categorical( p = ";
    if ( p != NULL )
    {
        if ( p->getName() != "" )
            o << p->getName();
        else
            o << "?";
    }
    else
    {
        o << "?";
    }
    o << " )";
}


/** Set a member variable */
void Dist_categorical::setConstParameter( const std::string& name, const RevPtr<const RevVariable>& var )
{
    if ( name == "p" )
    {
        p = var;
    }
    else
    {
        TypedDistribution<Natural>::setConstParameter( name, var );
    }
}
