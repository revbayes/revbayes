#include <cmath>
#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_BirthDeathJeffreys.h"
#include "MemberProcedure.h"
#include "BirthDeathJeffreysDistribution.h"
#include "ModelVector.h"
#include "Real.h"
#include "RealPos.h"
#include "RlMatrixRealSymmetric.h"
#include "ConstantNode.h"
#include "DagMemberFunction.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "MethodTable.h"
#include "ModelObject.h"
#include "Natural.h"
#include "OptionRule.h"
#include "RbException.h"
#include "RbHelpReference.h"
#include "RbVector.h"
#include "RevNullObject.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "RlDagMemberFunction.h"
#include "RlDeterministicNode.h"
#include "RlDistribution.h"
#include "RlStochasticNode.h"
#include "RlString.h"
#include "RlTypedDistribution.h"
#include "RlTypedFunction.h"
#include "RlUtils.h"
#include "StochasticNode.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

namespace RevBayesCore { class MatrixReal; }

using namespace RevLanguage;


/**
 * Default constructor.
 * 
 * The default constructor does nothing except allocating the object.
 */
Dist_BirthDeathJeffreys::Dist_BirthDeathJeffreys(void) : TypedDistribution<ModelVector<RealPos> >()
{

}


/**
 * Create a new internal distribution object.
 *
 * This function simply dynamically allocates a new internal distribution object that can be 
 * associated with the variable. The internal distribution object is created by calling its
 * constructor and passing the distribution-parameters (other DAG nodes) as arguments of the 
 * constructor. The distribution constructor takes care of the proper hook-ups.
 *
 * \return A new internal distribution object.
 */
RevBayesCore::BirthDeathJeffreysDistribution* Dist_BirthDeathJeffreys::createDistribution( void ) const
{

    // the start age
    RevBayesCore::TypedDagNode<double>* sa       = static_cast<const RealPos &>( start_age->getRevObject() ).getDagNode();

    // the start condition
    bool uo = ( start_condition == "originAge" ? true : false );

    // sampling condition
    const std::string& cond                     = static_cast<const RlString &>( condition->getRevObject() ).getValue();

    // sampling probability
    RevBayesCore::TypedDagNode<double>* r       = static_cast<const Probability &>( rho->getRevObject() ).getDagNode();

    // upper limit
    double lim = static_cast<const RealPos &>( limit->getRevObject() ).getValue();


    RevBayesCore::BirthDeathJeffreysDistribution* d = new RevBayesCore::BirthDeathJeffreysDistribution(sa, cond, uo, r, lim);

    return d;
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process. 
 */
Dist_BirthDeathJeffreys* Dist_BirthDeathJeffreys::clone( void ) const
{
    
    return new Dist_BirthDeathJeffreys(*this);
}

/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Dist_BirthDeathJeffreys::getClassType(void)
{ 
    
    static std::string rev_type = "Dist_BirthDeathJeffreys";
    
	return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Dist_BirthDeathJeffreys::getClassTypeSpec(void)
{ 
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution<ModelVector<Real> >::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_BirthDeathJeffreys::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "BirthDeathJeffreys";
    
    return d_name;
}


/** 
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the branch rate jump process are:
 * (1) the mean of the distribution.
 * (2) the standard deviation.
 *
 * \return The member rules.
 */
const MemberRules& Dist_BirthDeathJeffreys::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set ) 
    {

        std::vector<std::string> aliases;
        aliases.push_back("rootAge");
        aliases.push_back("originAge");
        dist_member_rules.push_back( new ArgumentRule( aliases, RealPos::getClassTypeSpec(), "The start time of the process.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    	dist_member_rules.push_back( new ArgumentRule( "rho",   Probability::getClassTypeSpec(), "The taxon sampling fraction.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Probability(1.0) ) );

        std::vector<std::string> optionsCondition;
        optionsCondition.push_back( "time" );
        optionsCondition.push_back( "survival" );
        optionsCondition.push_back( "nTaxa" );
        dist_member_rules.push_back( new OptionRule( "condition", new RlString("survival"), optionsCondition, "The condition of the process." ) );

    	dist_member_rules.push_back( new ArgumentRule( "limit", RealPos::getClassTypeSpec(), "The upper limit of speciation and extinction rates.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RealPos(10.0) ) );

        rules_set = true;
    }
    
    return dist_member_rules;

}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Dist_BirthDeathJeffreys::getTypeSpec( void ) const
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
void Dist_BirthDeathJeffreys::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{

    if ( name == "rootAge" || name == "originAge" )
	{
		start_age = var;
		start_condition = name;
	}
    else if ( name == "condition" )
    {
        condition = var;
    }
    else if ( name == "rho" )
    {
    	rho = var;
    }
    else if (name == "limit")
    {
    	limit = var;
    }
    else 
    {
        Distribution::setConstParameter(name, var);
    }
}
