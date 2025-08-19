
#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "EssMax.h"
#include "OptionRule.h"
#include "RlAbstractConvergenceStoppingRule.h"
#include "RbException.h"
#include "RlString.h"
#include "SemMin.h"
#include "TypeSpec.h"
#include "Natural.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlStoppingRule.h"

namespace RevBayesCore { class BurninEstimatorContinuous; }


using namespace RevLanguage;

/**
 * Default constructor.
 * Create the default instance.
 */
AbstractConvergenceStoppingRule::AbstractConvergenceStoppingRule(void) : StoppingRule()
{
    
}


RevBayesCore::BurninEstimatorContinuous* AbstractConvergenceStoppingRule::constructBurninEstimator( void )
{
    
    // create a new Burnin Estimator instance
    const std::string &bm = static_cast<const RlString &>( burninMethod->getRevObject() ).getValue();
    
    RevBayesCore::BurninEstimatorContinuous *burninEst = NULL;
    
    if ( bm == "ESS" )
    {
        burninEst = new RevBayesCore::EssMax();
    }
    else if ( bm == "SEM" )
    {
        burninEst = new RevBayesCore::SemMin();
    }
    else
    {
        throw RbException("Unknown burnin estimation method");
    }
    
    return burninEst;
}


/** Get Rev type of object */
const std::string& AbstractConvergenceStoppingRule::getClassType(void)
{
    
    static std::string rev_type = "AbstractConvergenceStoppingRule";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& AbstractConvergenceStoppingRule::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( StoppingRule::getClassTypeSpec() ) );
    
    return rev_type_spec;
}



/** Return member rules */
const MemberRules& AbstractConvergenceStoppingRule::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        memberRules.push_back( new ArgumentRule( "filename" , RlString::getClassTypeSpec(), "The name of the file containing the samples.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule( "frequency", Natural::getClassTypeSpec() , "The frequency how often to check for convergence.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(10000) ) );
        
        std::vector<std::string> bMethods;
        bMethods.push_back( "ESS" );
        bMethods.push_back( "SEM" );
        //        optionsUnits.push_back( "fixed" );
        memberRules.push_back( new OptionRule( "burninMethod", new RlString("ESS"), bMethods, "Which type of burnin method to use." ) );
        
        
        rules_set = true;
    }
    
    return memberRules;
}

/** Get type spec */
const TypeSpec& AbstractConvergenceStoppingRule::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Set a member variable */
void AbstractConvergenceStoppingRule::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "burninMethod" )
    {
        burninMethod = var;
    }
    else if ( name == "filename" )
    {
        filename = var;
    }
    else if ( name == "frequency" )
    {
        frequency = var;
    }
    else
    {
        StoppingRule::setConstParameter(name, var);
    }
    
}
