#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_PhylodynamicBDP.h"
#include "ModelVector.h"
#include "OptionRule.h"


namespace RevBayesCore { class Tree; }

using namespace RevLanguage;

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Dist_PhylodynamicBDP* Dist_PhylodynamicBDP::clone( void ) const
{
    return new Dist_PhylodynamicBDP(*this);
}

/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Dist_PhylodynamicBDP::getClassType( void )
{

    static std::string rev_type = "Dist_PhylodynamicBDSTP";

    return rev_type;
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_PhylodynamicBDP::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "PhylodynamicBDP" );
    a_names.push_back( "PhylodynamicBDSTP" );
    return a_names;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_PhylodynamicBDP::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "PhylodynamicBirthDeathProcess";

    return d_name;
}

/**
 * Get the member rules used to create the constructor of this object.
 *
 * Changes from the BDSTP default are:
 * (1) r is no longer required, defaults to 1 (but can be changed)
 * (2) rho/Phi is no longer required, defaults to 0 (but can be changed)
 * (3) No burst or mass extinctions
 *
 * \return The member rules.
 */
const MemberRules& Dist_PhylodynamicBDP::getParameterRules(void) const
{

    static MemberRules dist_member_rules;
    static bool rules_set = false;

    if ( rules_set == false )
    {
        Dist_BDSTP::addCommonRules(dist_member_rules);
        // no burst or mass extinctions

        std::vector<TypeSpec> event_sampling_paramTypes;
        event_sampling_paramTypes.push_back(Probability::getClassTypeSpec());
        event_sampling_paramTypes.push_back(ModelVector<Probability>::getClassTypeSpec());

        std::vector<std::string> aliases_event_sampling;
        aliases_event_sampling.push_back("Phi");
        aliases_event_sampling.push_back("rho");
        dist_member_rules.push_back(new ArgumentRule(aliases_event_sampling, event_sampling_paramTypes, "The probability of sampling taxa at sampling events (at present only if input is scalar).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Probability(0.0)));

        std::vector<TypeSpec> other_event_paramTypes;
        other_event_paramTypes.push_back(ModelVector<Probability>::getClassTypeSpec());
        dist_member_rules.push_back(new ArgumentRule("R", other_event_paramTypes, "The treatment probabilities for the sampling events (excluding sampling at present).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL));

        std::vector<TypeSpec> rTypes;
        rTypes.push_back(Probability::getClassTypeSpec());
        rTypes.push_back(ModelVector<Probability>::getClassTypeSpec());
        dist_member_rules.push_back(new ArgumentRule("r", rTypes, "The probabilit(y|ies) of death upon sampling (treatment).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Probability(1.0)));
        dist_member_rules.push_back(new ArgumentRule("rTimeline", ModelVector<RealPos>::getClassTypeSpec(), "The rate interval change times of the (serial) treatment probability.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL));

        rules_set = true;
    }

    return dist_member_rules;
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
void Dist_PhylodynamicBDP::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    if( name == "Lambda" || name == "Mu" || name == "LambdaTimeline" || name == "MuTimeline" ) {
        // do nothing - no bursts or mass extinctions in phylodynamics
    } else {
        Dist_BDSTP::setConstParameter(name, var);
    }
}
