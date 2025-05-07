#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "SoftBoundUniformNormalDistribution.h"
#include "Dist_softBoundUniformNormal.h"
#include "Probability.h"
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
Dist_SoftBoundUniformNormal::Dist_SoftBoundUniformNormal() : ContinuousDistribution()
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
RevBayesCore::SoftBoundUniformNormalDistribution* Dist_SoftBoundUniformNormal::createDistribution( void ) const
{
    
    // get the parameters
    RevBayesCore::TypedDagNode<double>* mi          = static_cast<const Real &>( min->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* ma          = static_cast<const Real &>( max->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* s           = static_cast<const RealPos &>( sd->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* p           = static_cast<const Probability &>( prob->getRevObject() ).getDagNode();
    RevBayesCore::SoftBoundUniformNormalDistribution*   d    = new RevBayesCore::SoftBoundUniformNormalDistribution(mi, ma, s, p);
    
    return d;
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Dist_SoftBoundUniformNormal* Dist_SoftBoundUniformNormal::clone( void ) const
{
    
    return new Dist_SoftBoundUniformNormal(*this);
}


/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Dist_SoftBoundUniformNormal::getClassType(void)
{
    
    static std::string rev_type = "Dist_SoftBoundUniformNormal";
    
    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Dist_SoftBoundUniformNormal::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( ContinuousDistribution::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_SoftBoundUniformNormal::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "SoftBoundUniformNormal";
    
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
const MemberRules& Dist_SoftBoundUniformNormal::getParameterRules(void) const
{
    
    static MemberRules distNormMemberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        distNormMemberRules.push_back( new ArgumentRule( "min", Real::getClassTypeSpec()       , "The min value of the uniform distribution.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY  ) );
        distNormMemberRules.push_back( new ArgumentRule( "max", Real::getClassTypeSpec()       , "The max value of the uniform distribution.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY  ) );
        distNormMemberRules.push_back( new ArgumentRule( "sd" , RealPos::getClassTypeSpec()    , "The standard deviation of the normal distribution.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY  ) );
        distNormMemberRules.push_back( new ArgumentRule( "p"  , Probability::getClassTypeSpec(), "The probability of being within the uniform distribution.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return distNormMemberRules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Dist_SoftBoundUniformNormal::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Print value for user */
void Dist_SoftBoundUniformNormal::printValue(std::ostream& o) const {
    
    o << "unifSBN(min=";
    if ( min != NULL ) {
        o << min->getName();
    }
    else
    {
        o << "?";
    }
    o << ", max=";
    if ( max != NULL )
    {
        o << max->getName();
    }
    else
    {
        o << "?";
    }
    o << ", sd=";
    if ( sd != NULL )
    {
        o << sd->getName();
    }
    else
    {
        o << "?";
    }
    o << ", p=";
    if ( prob != NULL )
    {
        o << prob->getName();
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
void Dist_SoftBoundUniformNormal::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "min" )
    {
        min = var;
    }
    else if ( name == "max" )
    {
        max = var;
    }
    else if ( name == "sd" )
    {
        sd = var;
    }
    else if ( name == "p" )
    {
        prob = var;
    }
    else
    {
        Distribution::setConstParameter(name, var);
    }
}
