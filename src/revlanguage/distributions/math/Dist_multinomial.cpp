#include <math.h>
#include <cstddef>
#include <cstdint>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_multinomial.h"
#include "MultinomialDistribution.h"
#include "RbVector.h"
#include "RlSimplex.h"
#include "ConstantNode.h"
#include "DagMemberFunction.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "ModelVector.h"
#include "Natural.h"
#include "RbHelpReference.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "RlDagMemberFunction.h"
#include "RlDeterministicNode.h"
#include "RlStochasticNode.h"
#include "RlTypedDistribution.h"
#include "RlTypedFunction.h"
#include "StochasticNode.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

namespace RevBayesCore { class Simplex; }

using namespace RevLanguage;

Dist_multinomial::Dist_multinomial(void) : TypedDistribution< ModelVector<Natural> >()
{
    
}


Dist_multinomial::~Dist_multinomial(void)
{
    
}



Dist_multinomial* Dist_multinomial::clone( void ) const
{
    
    return new Dist_multinomial(*this);
}


RevBayesCore::MultinomialDistribution* Dist_multinomial::createDistribution( void ) const
{
    
    // get the parameters
    RevBayesCore::TypedDagNode<std::int64_t>*                              vn = static_cast<const Natural     &>( n->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::Simplex >*          vp = static_cast<const Simplex &>( p->getRevObject() ).getDagNode();
    RevBayesCore::MultinomialDistribution*                        d  = new RevBayesCore::MultinomialDistribution( vp, vn );
    return d;
}



/* Get Rev type of object */
const std::string& Dist_multinomial::getClassType(void)
{
    
    static std::string rev_type = "Dist_multinomial";
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_multinomial::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution< ModelVector<Natural> >::getClassTypeSpec() ) );
    return rev_type_spec;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_multinomial::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "multinomial";
    
    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Dist_multinomial::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        dist_member_rules.push_back( new ArgumentRule( "p", Simplex::getClassTypeSpec(), "The simplex of probabilities for the categories.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "size", Natural::getClassTypeSpec(), "The number of draws.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    
    return dist_member_rules;
}


const TypeSpec& Dist_multinomial::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    return ts;
}


/** Print value for user */
void Dist_multinomial::printValue(std::ostream& o) const
{
    
    o << "multinomial(p=";
    if ( p != NULL )
        o << p->getName();
    else
        o << "?";
    o << ", size=";
    if ( n != NULL )
        o << n->getName();
    else
        o << "?";
    o << ")";
}


/** Set a member variable */
void Dist_multinomial::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "p" )
    {
        p = var;
    }
    else if ( name == "size" )
    {
        n = var;
    }
    else
    {
        TypedDistribution< ModelVector<Natural> >::setConstParameter(name, var);
    }
}
