#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_FBDP.h"
#include "ModelVector.h"
#include "OptionRule.h"
#include "RlConstantNode.h"

namespace RevBayesCore { class Tree; }

using namespace RevLanguage;

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Dist_FBDP* Dist_FBDP::clone( void ) const
{
    return new Dist_FBDP(*this);
}


/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Dist_FBDP::getClassType( void )
{

    static std::string rev_type = "Dist_FBDP";

    return rev_type;
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_FBDP::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "FBDP" );

    return a_names;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_FBDP::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "FossilizedBirthDeath";

    return d_name;
}

/**
 * Gets a fixed removal probability r = 0.
 */
RevBayesCore::DagNode* Dist_FBDP::getRemovalProbability( void ) const {
    return new ConstantNode<double>("r", new double(0.0) );
}

/**
 * Get the member rules used to create the constructor of this object.
 *
 * Changes from the BDSTP default are:
 * (1) r is no longer required, fixed to 0
 * (2) No mass sampling events
 *
 * \return The member rules.
 */
const MemberRules& Dist_FBDP::getParameterRules(void) const
{

    static MemberRules dist_member_rules;
    static bool rules_set = false;

    if ( rules_set == false )
    {
        Dist_BDSTP::addCommonRules(dist_member_rules);
        Dist_BDSTP::addBurstRules(dist_member_rules);

        // no removal probability - fixed to 0 for FBD
        std::vector<TypeSpec> event_sampling_paramTypes;
        event_sampling_paramTypes.push_back(Probability::getClassTypeSpec());
        event_sampling_paramTypes.push_back(ModelVector<Probability>::getClassTypeSpec());

        std::vector<std::string> aliases_event_sampling;
        aliases_event_sampling.push_back("rho");
        aliases_event_sampling.push_back("Phi");
        dist_member_rules.push_back(new ArgumentRule(aliases_event_sampling, event_sampling_paramTypes, "The probability of sampling taxa at sampling events (at present only if input is scalar).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY));

        rules_set = true;
    }

    return dist_member_rules;
}

/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Dist_FBDP::getTypeSpec( void ) const
{

    static TypeSpec ts = getClassTypeSpec();

    return ts;
}


/**
 * Set a member variable.
 *
 * Sets a member variable with the given name and store the pointer to the variable.
 * The value of the variable might still change but this function needs to be called again if the pointer to
 * the variable changes. The current values will be used to create the distribution object.
 *
 * \param[in]    name     Name of the member variable.
 * \param[in]    var      Pointer to the variable.
 */
void Dist_FBDP::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    if ( name == "r" || name == "rTimeline" || name == "R" ) {
        // do nothing - no removal in the FBD
    } else {
        Dist_BDSTP::setConstParameter(name, var);
    }

}
