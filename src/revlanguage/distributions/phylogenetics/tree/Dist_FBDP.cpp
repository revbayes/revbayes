#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_FBDP.h"
#include "ModelVector.h"
#include "OptionRule.h"
#include "BirthDeathSamplingTreatmentProcess.h"
#include "Probability.h"
#include "RealPos.h"
#include "RlString.h"
#include "RlTaxon.h"
#include "RlTimeTree.h"
#include "AbstractBirthDeathProcess.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "RbVector.h"
#include "RevNullObject.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlBirthDeathProcess.h"
#include "RlConstantNode.h"
#include "Taxon.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

namespace RevBayesCore { class Tree; }

using namespace RevLanguage;

/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Dist_FBDP::Dist_FBDP() : Dist_BDSTP()
{   
    const RevPtr t_s = new RevVariable(new Real(0.0) , "r");
    Dist_BDSTP::setParameter("r", t_s);
}


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


void Dist_FBDP::addSamplingAndRemovalRules(MemberRules& dist_member_rules) const
{    
    std::vector<TypeSpec> event_sampling_paramTypes;
    event_sampling_paramTypes.push_back( Probability::getClassTypeSpec() );
    event_sampling_paramTypes.push_back( ModelVector<Probability>::getClassTypeSpec() );
    
    std::vector<std::string> aliases_event_sampling;
    aliases_event_sampling.push_back("rho");
    dist_member_rules.push_back( new ArgumentRule( aliases_event_sampling,     event_sampling_paramTypes, "The probability of sampling taxa at sampling events (at present only if input is scalar).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
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
