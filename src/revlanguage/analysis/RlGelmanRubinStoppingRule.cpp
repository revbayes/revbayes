
#include <cstddef>
#include <ostream>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "GelmanRubinStoppingRule.h"
#include "RlGelmanRubinStoppingRule.h"
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
GelmanRubinStoppingRule::GelmanRubinStoppingRule(void) : AbstractConvergenceStoppingRule()
{
    
}


/**
 * Clone object
 */
GelmanRubinStoppingRule* GelmanRubinStoppingRule::clone(void) const
{
    
    return new GelmanRubinStoppingRule(*this);
}


void GelmanRubinStoppingRule::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new stopping rule
    double r = static_cast<const RealPos &>( R->getRevObject() ).getValue();
    int fq = (int)static_cast<const Natural &>( frequency->getRevObject() ).getValue();
    const std::string &fn = static_cast<const RlString &>( filename->getRevObject() ).getValue();
    
    RevBayesCore::BurninEstimatorContinuous *burninEst = constructBurninEstimator();
    
    value = new RevBayesCore::GelmanRubinStoppingRule(r, fn, size_t(fq), burninEst);
}


/** Get Rev type of object */
const std::string& GelmanRubinStoppingRule::getClassType(void)
{
    
    static std::string rev_type = "GelmanRubinStoppingRule";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& GelmanRubinStoppingRule::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( AbstractConvergenceStoppingRule::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string GelmanRubinStoppingRule::getConstructorFunctionName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "srGelmanRubin";
    
    return c_name;
}


/** Return member rules */
const MemberRules& GelmanRubinStoppingRule::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        memberRules.push_back( new ArgumentRule( "R", RealPos::getClassTypeSpec(), "The maximum allowed potential scale reduction factor.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = AbstractConvergenceStoppingRule::getParameterRules();
        memberRules.insert( memberRules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return memberRules;
}

/** Get type spec */
const TypeSpec& GelmanRubinStoppingRule::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void GelmanRubinStoppingRule::printValue(std::ostream &o) const
{
    
    o << "GelmanRubinStoppingRule";
}


/** Set a member variable */
void GelmanRubinStoppingRule::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "R" )
    {
        R = var;
    }
    else
    {
        AbstractConvergenceStoppingRule::setConstParameter(name, var);
    }
    
}
