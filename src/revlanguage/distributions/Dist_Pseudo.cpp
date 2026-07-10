#include "Dist_Pseudo.h"

RevBayesCore::PseudoDist* Dist_Pseudo::createDistribution( void ) const
{
    auto pdl_node     = static_cast<const PseudoDataLikelihood&>( pseudoDataLikelihood->getRevObject() ).getDagNode();

    return new RevBayesCore::PseudoDist(pdl_node);
}

const MemberRules&   Dist_Pseudo::getParameterRules(void) const
{
    static MemberRules dist_member_rules;
    static bool rules_set = false;

    if ( rules_set == false )
    {
        dist_member_rules.push_back( new ArgumentRule( "pseudoDataLikelihood", PseudoDataLikelihood::getClassTypeSpec(), "The likelihood of the unspecified data under the unspecified distribution.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return dist_member_rules;
}


void Dist_Pseudo::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
        if ( name == "pseudoDataLikelihood" )
        { 
            pseudoDataLikelihood = var;
        }
        else
        {
            TypedDistribution< PseudoObservation >::setConstParameter(name, var);
        }
}



