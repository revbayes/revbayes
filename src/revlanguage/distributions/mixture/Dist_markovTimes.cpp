#include <cmath>
#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_markovTimes.h"

#include "OrderedEventTimes.h"
#include "ModelVector.h"
#include "Natural.h"
#include "RlDistributionMemberFunction.h"
#include "RlString.h"
#include "StochasticNode.h"
#include "WorkspaceVector.h"
#include "ConstantNode.h"
#include "DagMemberFunction.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DistributionMemberFunction.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "MethodTable.h"
#include "ModelObject.h"
#include "MultiValueEvent.h"
#include "RbHelpReference.h"
#include "RbVector.h"
#include "Real.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "RlDagMemberFunction.h"
#include "RlDeterministicNode.h"
#include "RlDistribution.h"
#include "RlOrderedEventTimes.h"
#include "RlStochasticNode.h"
#include "RlTypedDistribution.h"
#include "RlTypedFunction.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"
#include "WorkspaceToCoreWrapperObject.h"

using namespace RevLanguage;


/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Dist_markovTimes::Dist_markovTimes() : TypedDistribution< RlOrderedEventTimes >()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
Dist_markovTimes* Dist_markovTimes::clone( void ) const
{
    
    return new Dist_markovTimes(*this);
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
RevBayesCore::MarkovTimesDistribution* Dist_markovTimes::createDistribution( void ) const
{
	// get the parameters
    RevBayesCore::TypedDagNode<double>* the_rate = static_cast<const RealPos &>( rate->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* the_age  = static_cast<const RealPos &>( age->getRevObject() ).getDagNode();

    // make the distribution
	RevBayesCore::MarkovTimesDistribution* d = new RevBayesCore::MarkovTimesDistribution(the_rate, the_age);

	return d;
}



/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Dist_markovTimes::getClassType(void)
{
    
    static std::string rev_type = "Dist_markovTimes";
    
    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Dist_markovTimes::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution<RlOrderedEventTimes>::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_markovTimes::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    
    return a_names;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_markovTimes::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "MarkovTimes";
    
    return d_name;
}

MethodTable Dist_markovTimes::getDistributionMethods( void ) const
{
    MethodTable methods = TypedDistribution<RlOrderedEventTimes>::getDistributionMethods();

    return methods;
}

/**
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the MultiValueEvent distribution are:
 * (1) the rate lambda which must be a positive real between 0 and 1 (= a probability).
 *
 * \return The member rules.
 */
const MemberRules& Dist_markovTimes::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        
        dist_member_rules.push_back( new ArgumentRule( "rate", RealPos::getClassTypeSpec(), "The rate of the process.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "age",  RealPos::getClassTypeSpec(), "The age of the process.",  ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }
    
    return dist_member_rules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Dist_markovTimes::getTypeSpec( void ) const
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
void Dist_markovTimes::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
     
    if ( name == "rate" )
    {
        rate = var;
    }
    else if ( name == "age" )
    {
        age = var;
    }
    else
    {
        TypedDistribution< RlOrderedEventTimes >::setConstParameter(name, var);
    }
    
}
