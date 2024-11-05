#include "Func_time.h"

#include <string>

#include "Natural.h"
#include "OptionRule.h"
#include "RevNullObject.h"
#include "TypeSpec.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "RbHelpReference.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "RlString.h"

#include "boost/date_time/posix_time/posix_time.hpp" // IWYU pragma: keep
#include "boost/date_time/gregorian/gregorian.hpp" // IWYU pragma: keep
#include "boost/date_time/gregorian/greg_month.hpp" // IWYU pragma: keep

using namespace RevLanguage;

Func_time::Func_time() : Procedure()
{
    
}

/* Clone object */
Func_time* Func_time::clone( void ) const
{
    
    return new Func_time( *this );
}


/** Execute function: We rely on getValue and overloaded push_back to provide functionality */
RevPtr<RevVariable> Func_time::execute( void )
{
    
    const std::string &option = static_cast<const RlString &>( args[0].getVariable()->getRevObject() ).getValue();
    
    boost::posix_time::ptime t = boost::posix_time::microsec_clock::local_time();

    long std::int64_t = 0;
    if ( option == "year" )
    {
        time = t.date().year();
    }
    else if ( option == "day" )
    {
        time = t.date().day_of_year();
    }
    else if ( option == "seconds" )
    {
        time = t.time_of_day().total_seconds();
    }
    else if ( option == "milliseconds" )
    {
        time = t.time_of_day().total_milliseconds();
    }
    else if ( option == "fromBeginning" )
    {
        boost::posix_time::ptime epoch_start(boost::posix_time::ptime(boost::gregorian::date(1970,1,1)));
        boost::posix_time::time_duration duration = t - epoch_start;
         
        time = duration.total_milliseconds();
    }
    
    return RevPtr<RevVariable>( new RevVariable( new Natural( time ) ) );
}


/** Get argument rules */
const ArgumentRules& Func_time::getArgumentRules( void ) const
{
    
    static ArgumentRules argument_rules = ArgumentRules();
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        std::vector<std::string> options;
        options.push_back( "year" );
        options.push_back( "day" );
        options.push_back( "seconds" );
        options.push_back( "milliseconds" );
        options.push_back( "fromBeginning" );
        argument_rules.push_back( new OptionRule( "option", new RlString("fromBeginning"), options, "The format of the time." ) );
        
        rules_set = true;
    }
    
    return argument_rules;
}


/** Get Rev type of object */
const std::string& Func_time::getClassType(void)
{
    
    static std::string rev_type = "Func_time";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_time::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_time::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "time";
    
    return f_name;
}



/** Get type spec */
const TypeSpec& Func_time::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_time::getReturnType( void ) const
{
    
    return RevNullObject::getClassTypeSpec();
}
