
#include <cstddef>
#include <ostream>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MinEssStoppingRule.h"
#include "RlMinEssStoppingRule.h"
#include "RealPos.h"
#include "RlString.h"
#include "TypeSpec.h"
#include "Natural.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlAbstractConvergenceStoppingRule.h"
#include "StoppingRule.h"

namespace RevBayesCore { class BurninEstimatorContinuous; }


using namespace RevLanguage;

/**
 * Default constructor.
 * Create the default instance.
 */
MinEssStoppingRule::MinEssStoppingRule(void) : AbstractConvergenceStoppingRule()
{
    
}


/**
 * Clone object
 */
MinEssStoppingRule* MinEssStoppingRule::clone(void) const
{
    
    return new MinEssStoppingRule(*this);
}


void MinEssStoppingRule::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new stopping rule
    double min = static_cast<const RealPos &>( minEss->getRevObject() ).getValue();
    int fq = static_cast<const Natural &>( frequency->getRevObject() ).getValue();
    const std::string &fn = static_cast<const RlString &>( filename->getRevObject() ).getValue();
    
    RevBayesCore::BurninEstimatorContinuous *burninEst = constructBurninEstimator();
    
    value = new RevBayesCore::MinEssStoppingRule(min, fn, size_t(fq), burninEst);

}


/** Get Rev type of object */
const std::string& MinEssStoppingRule::getClassType(void)
{
    
    static std::string rev_type = "MinEssStoppingRule";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& MinEssStoppingRule::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( AbstractConvergenceStoppingRule::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string MinEssStoppingRule::getConstructorFunctionName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "srMinESS";
    
    return c_name;
}


/** Return member rules */
const MemberRules& MinEssStoppingRule::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        memberRules.push_back( new ArgumentRule( "minEss"   , RealPos::getClassTypeSpec(), "The minimum ESS threshold when stopping is allowed.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = AbstractConvergenceStoppingRule::getParameterRules();
        memberRules.insert( memberRules.end(), inheritedRules.begin(), inheritedRules.end() );


        rules_set = true;
    }
    
    return memberRules;
}

/** Get type spec */
const TypeSpec& MinEssStoppingRule::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void MinEssStoppingRule::printValue(std::ostream &o) const
{
    
    o << "MinEssStoppingRule";
}


/** Set a member variable */
void MinEssStoppingRule::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "minEss" )
    {
        minEss = var;
    }
    else
    {
        AbstractConvergenceStoppingRule::setConstParameter(name, var);
    }
    
}
