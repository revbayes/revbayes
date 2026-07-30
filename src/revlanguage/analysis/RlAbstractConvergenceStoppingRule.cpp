
#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "EssMax.h"
#include "FixedBurnin.h"
#include "OptionRule.h"
#include "Probability.h"
#include "RlAbstractConvergenceStoppingRule.h"
#include "RbException.h"
#include "RlString.h"
#include "RlUserInterface.h"
#include "SemMin.h"
#include "TypeSpec.h"
#include "IntegerPos.h"
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
        // We want to throw a warning when the user explicitly specifies the `burnin` argument while simultaneously setting
        // `burninMethod` to "ESS" or "SEM". However, the argument has a default value, so it will not be NULL even if the user does
        // not explicitly set it. To circumvent the issue, we will exploit an implicit convention in Function::processArguments():
        // when an optional parameter is assigned its ArgumentRule default, its name is prepended with a leading dot.
        if ( burnin != NULL && burnin->getName() != ".burnin" )
        {
            RBOUT( "Warning: the `burnin` argument is ignored when `burninMethod` is \"ESS\"." );
        }
        burninEst = new RevBayesCore::EssMax();
    }
    else if ( bm == "SEM" )
    {
        if ( burnin != NULL && burnin->getName() != ".burnin" )
        {
            RBOUT( "Warning: the `burnin` argument is ignored when `burninMethod` is \"SEM\"." );
        }
        burninEst = new RevBayesCore::SemMin();
    }
    else if ( bm == "fixed" )
    {
        double fraction = static_cast<const Probability &>( burnin->getRevObject() ).getValue();
        burninEst = new RevBayesCore::FixedBurnin( fraction );
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
        memberRules.push_back( new ArgumentRule( "frequency", IntegerPos::getClassTypeSpec() , "The frequency how often to check for convergence.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new IntegerPos(10000) ) );
        
        std::vector<std::string> bMethods;
        bMethods.push_back( "ESS" );
        bMethods.push_back( "SEM" );
        bMethods.push_back( "fixed" );
        memberRules.push_back( new OptionRule( "burninMethod", new RlString("ESS"), bMethods, "Which type of burnin method to use." ) );
        memberRules.push_back( new ArgumentRule( "burnin", Probability::getClassTypeSpec(), "The fraction of samples to discard as burnin (only used when burninMethod is \"fixed\").", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Probability(0.25) ) );
        
        
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
    
    if ( name == "burnin" )
    {
        burnin = var;
    }
    else if ( name == "burninMethod" )
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
