#include "Dist_Require.h"


#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "GenericFunction.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RbException.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "RlBoolean.h"
#include "TypeSpec.h"

namespace Core = RevBayesCore;

RevBayesCore::RequireDist* Dist_Require::createDistribution( void ) const
{
    // get the parameters
    auto pred  = static_cast<const RlBoolean&>( predicate->getRevObject() ).getDagNode();
    auto w     = static_cast<const RealPos&>( weight->getRevObject() ).getDagNode();

    return new RevBayesCore::RequireDist(pred, w);
}

const MemberRules&   Dist_Require::getParameterRules(void) const
{
    static MemberRules dist_member_rules;
    static bool rules_set = false;

    if ( rules_set == false )
    {
        dist_member_rules.push_back( new ArgumentRule( "predicate", RlBoolean::getClassTypeSpec(), "The predicate.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "weight", RealPos::getClassTypeSpec(), "The odds ratio in favor of the predicate.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(std::numeric_limits<double>::infinity() ) ) );

        rules_set = true;
    }

    return dist_member_rules;
}


void Dist_Require::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    if ( name == "predicate" )
    { 
        predicate = var;
    }
    else if ( name == "weight" )
    { 
        weight = var;
    }
    else
    {
        TypedDistribution< PseudoObservation >::setConstParameter(name, var);
    }
}

