#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "OrnsteinUhlenbeckProcess.h"
#include "Process_OrnsteinUhlenbeck.h"
#include "Real.h"
#include "RealPos.h"
#include "RbHelpReference.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlContinuousDistribution.h"
#include "RlDistribution.h"
#include "TypeSpec.h"

namespace RevBayesCore { template <class valueType> class TypedDagNode; }

using namespace RevLanguage;


/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
OrnsteinUhlenbeckProcess::OrnsteinUhlenbeckProcess() : ContinuousDistribution()
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
RevBayesCore::OrnsteinUhlenbeckProcess* OrnsteinUhlenbeckProcess::createDistribution( void ) const
{
    
    // get the parameters
    RevBayesCore::TypedDagNode<double>* x           = static_cast<const Real &>( x0->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* th          = static_cast<const Real &>( theta->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* s           = static_cast<const RealPos &>( sigma->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* a           = static_cast<const RealPos &>( alpha->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* ti          = static_cast<const RealPos &>( time->getRevObject() ).getDagNode();
    RevBayesCore::OrnsteinUhlenbeckProcess*   d     = new RevBayesCore::OrnsteinUhlenbeckProcess(x, th, a, s, ti);
    
    return d;
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
OrnsteinUhlenbeckProcess* OrnsteinUhlenbeckProcess::clone( void ) const
{
    
    return new OrnsteinUhlenbeckProcess(*this);
}


/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& OrnsteinUhlenbeckProcess::getClassType(void)
{
    
    static std::string rev_type = "OrnsteinUhlenbeckProcess";
    
    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& OrnsteinUhlenbeckProcess::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( ContinuousDistribution::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> OrnsteinUhlenbeckProcess::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "OU" );
    
    return a_names;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string OrnsteinUhlenbeckProcess::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "OrnsteinUhlenbeck";
    
    return d_name;
}


/**
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the branch rate jump process are:
 * (1) the min of the uniform distribution.
 * (2) the max of the uniform distribution.
 * (3) the standard deviation of the normal.
 * (4) the probability of belonging to distribution one.
 *
 * \return The member rules.
 */
const MemberRules& OrnsteinUhlenbeckProcess::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        dist_member_rules.push_back( new ArgumentRule( "x0"   , Real::getClassTypeSpec()   , "The root parameter value.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY  ) );
        dist_member_rules.push_back( new ArgumentRule( "theta", Real::getClassTypeSpec()   , "The location of the optimum parameter.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY  ) );
        dist_member_rules.push_back( new ArgumentRule( "alpha", RealPos::getClassTypeSpec(), "The attraction to the optimum parameter.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY  ) );
        dist_member_rules.push_back( new ArgumentRule( "sigma", RealPos::getClassTypeSpec(), "The scaling parameter of the time.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "time" , RealPos::getClassTypeSpec(), "The duration of the process.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return dist_member_rules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& OrnsteinUhlenbeckProcess::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Print value for user */
void OrnsteinUhlenbeckProcess::printValue(std::ostream& o) const
{
    
    o << "OU(x0=";
    if ( x0 != NULL ) {
        o << x0->getName();
    }
    else
    {
        o << "?";
    }
    o << ", theta=";
    if ( theta != NULL )
    {
        o << theta->getName();
    }
    else
    {
        o << "?";
    }
    o << ", alpha=";
    if ( alpha != NULL )
    {
        o << alpha->getName();
    }
    else
    {
        o << "?";
    }
    o << ", sigma=";
    if ( sigma != NULL )
    {
        o << sigma->getName();
    }
    else
    {
        o << "?";
    }
    o << ", time=";
    if ( time != NULL )
    {
        o << time->getName();
    }
    else
    {
        o << "?";
    }
    o << ")";
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
void OrnsteinUhlenbeckProcess::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "x0" )
    {
        x0 = var;
    }
    else if ( name == "theta" )
    {
        theta = var;
    }
    else if ( name == "alpha" )
    {
        alpha = var;
    }
    else if ( name == "sigma" )
    {
        sigma = var;
    }
    else if ( name == "time" )
    {
        time = var;
    }
    else
    {
        Distribution::setConstParameter(name, var);
    }
}
