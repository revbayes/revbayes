#include "Func_profile.h"

#include <string>

#include "Real.h"
#include "OptionRule.h"
#include "RevNullObject.h"
#include "TypeSpec.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "RbHelpReference.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "RlString.h"

#include "UniqueProfiler.h"

using namespace RevLanguage;

Func_profile::Func_profile() : Procedure()
{

}

/* Clone object */
Func_profile* Func_profile::clone( void ) const
{

    return new Func_profile( *this );
}


/** Execute function: We rely on getValue and overloaded push_back to provide functionality */
RevPtr<RevVariable> Func_profile::execute( void )
{

    const std::string &command = static_cast<const RlString &>( args[0].getVariable()->getRevObject() ).getValue();
    const std::string &event = static_cast<const RlString &>( args[0].getVariable()->getRevObject() ).getValue();

    double duration = 0.;
    if ( command == "start" )
    {
        Utils::Profiling::UniqueProfiler::getInstance()->startEvent(event);
    }
    else if ( command == "stop" )
    {
        Utils::Profiling::UniqueProfiler::getInstance()->endEvent(event);
    }
    else if ( command == "reportTotalTime" )
    {
        duration = Utils::Profiling::UniqueProfiler::getInstance()->reportEventTotalDuration(event);
    }
    else if ( command == "reportPerCallAvgTime" )
    {
        duration = Utils::Profiling::UniqueProfiler::getInstance()->reportEventAvgPerCallDuration(event);
    }

    return new RevVariable( new Real( duration ) );
}


/** Get argument rules */
const ArgumentRules& Func_profile::getArgumentRules( void ) const
{

    static ArgumentRules argument_rules = ArgumentRules();
    static bool rules_set = false;

    if ( !rules_set )
    {
        std::vector<std::string> options;
        options.push_back( "start" );
        options.push_back( "stop" );
        options.push_back( "reportTotalTime" );
        options.push_back( "reportPerCallAvgTime" );

        argument_rules.push_back( new OptionRule( "command", new RlString("start"), options, "The command to apply (start, stop, reportTotalTime or reportPerCallAvgTime)." ) );
        argument_rules.push_back( new ArgumentRule( "event", RlString::getClassTypeSpec(), "The name of the event to track.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return argument_rules;
}


/** Get Rev type of object */
const std::string& Func_profile::getClassType(void)
{

    static std::string rev_type = "Func_profile";

    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_profile::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_profile::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "profile";

    return f_name;
}



/** Get type spec */
const TypeSpec& Func_profile::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}


/** Get return type */
const TypeSpec& Func_profile::getReturnType( void ) const
{

    return RevNullObject::getClassTypeSpec();
}
