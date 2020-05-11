#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ContinuousStochasticNode.h"
#include "Dist_pointMassPositive.h"
#include "Real.h"
#include "RealPos.h"
#include "PointMassDistribution.h"

using namespace RevLanguage;

Dist_pointMassPositive::Dist_pointMassPositive() : PositiveContinuousDistribution()
{

}


Dist_pointMassPositive::~Dist_pointMassPositive() {

}



Dist_pointMassPositive* Dist_pointMassPositive::clone( void ) const
{
    return new Dist_pointMassPositive (*this);
}


RevBayesCore::PointMassDistribution* Dist_pointMassPositive::createDistribution( void ) const
{
    // get the parameters
    RevBayesCore::TypedDagNode<double>* l   = static_cast<const RealPos &>( value->getRevObject() ).getDagNode();
    RevBayesCore::PointMassDistribution* d    = new RevBayesCore::PointMassDistribution(l);

    return d;
}



/* Get Rev type of object */
const std::string& Dist_pointMassPositive::getClassType(void)
{

    static std::string rev_type = "Dist_pointMass";

	return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_pointMassPositive::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( PositiveContinuousDistribution::getClassTypeSpec() ) );

	return rev_type_spec;
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_pointMassPositive::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "DiracDelta" );

    return a_names;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_pointMassPositive::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "PointMass";

    return d_name;
}


/** Return member rules */
const MemberRules& Dist_pointMassPositive::getParameterRules(void) const
{

    static MemberRules distUnifMemberRules;
    static bool rules_set = false;

    if ( !rules_set )
    {

        distUnifMemberRules.push_back( new ArgumentRule( "value", RealPos::getClassTypeSpec(), "The value where mass is concentrated.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return distUnifMemberRules;
}


const TypeSpec& Dist_pointMassPositive::getTypeSpec( void ) const
{

    static TypeSpec ts = getClassTypeSpec();

    return ts;
}


/** Print value for user */
void Dist_pointMassPositive::printValue(std::ostream& o) const
{

    o << " pointMass (value=";
    if ( value != NULL ) {
        o << value->getName();
    } else {
        o << "?";
    }
    o << ")";
}


/** Set a member variable */
void Dist_pointMassPositive::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{

    if ( name == "value" )
    {
        value = var;
    }
    else
    {
        PositiveContinuousDistribution::setConstParameter(name, var);
    }
}
