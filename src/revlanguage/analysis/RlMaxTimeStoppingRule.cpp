
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MaxTimeStoppingRule.h"
#include "RlMaxTimeStoppingRule.h"
#include "RealPos.h"
#include "RbException.h"
#include "OptionRule.h"
#include "RlString.h"
#include "TypeSpec.h"


using namespace RevLanguage;

/**
 * Default constructor.
 * Create the default instance.
 */
MaxTimeStoppingRule::MaxTimeStoppingRule(void) : StoppingRule()
{
    
}


/** 
 * Clone object 
 */
MaxTimeStoppingRule* MaxTimeStoppingRule::clone(void) const
{
    
    return new MaxTimeStoppingRule(*this);
}


void MaxTimeStoppingRule::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new stopping rule
    double mt = static_cast<const RealPos &>( maxTime->getRevObject() ).getValue();
    const std::string &u = static_cast<const RlString &>( unit->getRevObject() ).getValue();
    if ( u == "minutes" )
    {
        mt *= 60;
    }
    else if ( u == "hours" )
    {
        mt *= 3600;
    }
    
    value = new RevBayesCore::MaxTimeStoppingRule(mt);
}


/** Get Rev type of object */
const std::string& MaxTimeStoppingRule::getClassType(void)
{
    
    static std::string rev_type = "MaxTimeStoppingRule";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& MaxTimeStoppingRule::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( StoppingRule::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string MaxTimeStoppingRule::getConstructorFunctionName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "srMaxTime";
    
    return c_name;
}


/** Return member rules */
const MemberRules& MaxTimeStoppingRule::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        memberRules.push_back( new ArgumentRule("maxTime" , RealPos::getClassTypeSpec(), "The maximum time to run.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        
        std::vector<std::string> optionsUnits;
        optionsUnits.push_back( "seconds" );
        optionsUnits.push_back( "minutes" );
        optionsUnits.push_back( "hours" );
        memberRules.push_back( new OptionRule( "unit", new RlString("seconds"), optionsUnits, "The unit in which we measure the maximum time." ) );

        rules_set = true;
    }
    
    return memberRules;
}

/** Get type spec */
const TypeSpec& MaxTimeStoppingRule::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void MaxTimeStoppingRule::printValue(std::ostream &o) const
{
    
    o << "MaxTimeStoppingRule";
}


/** Set a member variable */
void MaxTimeStoppingRule::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
    if ( name == "maxTime" )
    {
        maxTime = var;
    }
    else if ( name == "unit" )
    {
        unit = var;
    }
    else
    {
        StoppingRule::setConstParameter(name, var);
    }
    
}
