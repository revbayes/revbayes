
#include <cstddef>
#include <ostream>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "StationarityStoppingRule.h"
#include "RlStationarityStoppingRule.h"
#include "Probability.h"
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
StationarityStoppingRule::StationarityStoppingRule(void) : AbstractConvergenceStoppingRule()
{
    
}


/**
 * Clone object
 */
StationarityStoppingRule* StationarityStoppingRule::clone(void) const
{
    
    return new StationarityStoppingRule(*this);
}


void StationarityStoppingRule::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new stopping rule
    double p = static_cast<const Probability &>( prob->getRevObject() ).getValue();
    int fq = static_cast<const Natural &>( frequency->getRevObject() ).getValue();
    const std::string &fn = static_cast<const RlString &>( filename->getRevObject() ).getValue();
    
    RevBayesCore::BurninEstimatorContinuous *burninEst = constructBurninEstimator();
    
    value = new RevBayesCore::StationarityStoppingRule(p, fn, size_t(fq), burninEst);
}


/** Get Rev type of object */
const std::string& StationarityStoppingRule::getClassType(void)
{
    
    static std::string rev_type = "StationarityStoppingRule";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& StationarityStoppingRule::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( AbstractConvergenceStoppingRule::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string StationarityStoppingRule::getConstructorFunctionName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "srStationarity";
    
    return c_name;
}


/** Return member rules */
const MemberRules& StationarityStoppingRule::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        memberRules.push_back( new ArgumentRule( "prob", Probability::getClassTypeSpec() , "The significance level.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = AbstractConvergenceStoppingRule::getParameterRules();
        memberRules.insert( memberRules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return memberRules;
}

/** Get type spec */
const TypeSpec& StationarityStoppingRule::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void StationarityStoppingRule::printValue(std::ostream &o) const
{
    
    o << "StationarityStoppingRule";
}


/** Set a member variable */
void StationarityStoppingRule::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "prob" )
    {
        prob = var;
    }
    else
    {
        AbstractConvergenceStoppingRule::setConstParameter(name, var);
    }
    
}
