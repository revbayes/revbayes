#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "BimodalLognormalDistribution.h"
#include "Dist_bimodalLnorm.h"
#include "Probability.h"
#include "Real.h"
#include "RealPos.h"
#include "RbHelpReference.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlDistribution.h"
#include "RlPositiveContinuousDistribution.h"
#include "TypeSpec.h"

namespace RevBayesCore { template <class valueType> class TypedDagNode; }

using namespace RevLanguage;


/**
 * Default constructor.
 * 
 * The default constructor does nothing except allocating the object.
 */
Dist_bimodalLnorm::Dist_bimodalLnorm() : PositiveContinuousDistribution()
{
    
    setGuiDistributionName("Bimodal LogNormal");
    setGuiDistributionToolTip("");
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
RevBayesCore::BimodalLognormalDistribution* Dist_bimodalLnorm::createDistribution( void ) const
{
    
    // get the parameters
    RevBayesCore::TypedDagNode<double>* m1          = static_cast<const Real &>( mean1->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* m2          = static_cast<const Real &>( mean2->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* s1          = static_cast<const RealPos &>( sd1->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* s2          = static_cast<const RealPos &>( sd2->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* p           = static_cast<const Probability &>( prob->getRevObject() ).getDagNode();
    RevBayesCore::BimodalLognormalDistribution*   d = new RevBayesCore::BimodalLognormalDistribution(m1, m2, s1, s2, p);
    
    return d;
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process. 
 */
Dist_bimodalLnorm* Dist_bimodalLnorm::clone( void ) const 
{
    
    return new Dist_bimodalLnorm(*this);
}


/**
 * Get Rev type of object 
 *
 * \return The class' name.
 */
const std::string& Dist_bimodalLnorm::getClassType(void) 
{ 
    
    static std::string rev_type = "Dist_bimodalLnorm";
    
	return rev_type; 
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Dist_bimodalLnorm::getClassTypeSpec(void) 
{ 
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( PositiveContinuousDistribution::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_bimodalLnorm::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "BimodalLognormal";
    
    return d_name;
}


/** 
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the branch rate jump process are:
 * (1) the mean of the first distribution.
 * (2) the mean of the second distribution.
 * (3) the standard first deviation.
 * (4) the standard second deviation.
 * (5) the probability of belonging to distribution one.
 *
 * \return The member rules.
 */
const MemberRules& Dist_bimodalLnorm::getParameterRules(void) const 
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set ) 
    {
        
        memberRules.push_back( new ArgumentRule( "mean1", Real::getClassTypeSpec()       , "The mean (in log-space) of the first lognormal distribution.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule( "mean2", Real::getClassTypeSpec()       , "The mean (in log-space) of the second lognormal distribution.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule( "sd1"  , RealPos::getClassTypeSpec()    , "The standard deviation (in log-space) of the first lognormal distribution.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule( "sd2"  , RealPos::getClassTypeSpec()    , "The standard deviation (in log-space) of the second lognormal distribution.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule( "p"    , Probability::getClassTypeSpec(), "The probability to belong to the first distribution.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return memberRules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Dist_bimodalLnorm::getTypeSpec( void ) const 
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Print value for user */
void Dist_bimodalLnorm::printValue(std::ostream& o) const {
    
    o << "blnorm(mean1=";
    if ( mean1 != NULL ) {
        o << mean1->getName();
    } 
    else 
    {
        o << "?";
    }
    o << ", mean2=";
    if ( mean2 != NULL ) 
    {
        o << mean2->getName();
    } 
    else 
    {
        o << "?";
    }
    o << ", sd1=";
    if ( sd1 != NULL ) 
    {
        o << sd1->getName();
    } 
    else 
    {
        o << "?";
    }
    o << ", sd2=";
    if ( sd2 != NULL ) 
    {
        o << sd2->getName();
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
void Dist_bimodalLnorm::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) 
{
 
    if ( name == "mean1" ) 
    {
        mean1 = var;
    }    
    else if ( name == "mean2" ) 
    {
        mean2 = var;
    }
    else if ( name == "sd1" ) 
    {
        sd1 = var;
    }
    else if ( name == "sd2" ) 
    {
        sd2 = var;
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
