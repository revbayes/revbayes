#include "Func_time.h"
#include "ArgumentRule.h"
#include "DeterministicNode.h"
#include "Ellipsis.h"
#include "Natural.h"
#include "OptionRule.h"
#include "RandomNumberGenerator.h"
#include "RandomNumberFactory.h"
#include "RbUtil.h"
#include "RevNullObject.h"
#include "RealPos.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"

#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/date_time/gregorian/gregorian.hpp"
#include "boost/date_time/gregorian/greg_month.hpp"

#include <ctime>
#include <ostream>

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

    long time = 0;
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
        boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
        boost::posix_time::ptime t2(boost::posix_time::min_date_time);
        boost::posix_time::time_duration duration = t1-t2;
         
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


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Func_time::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Sebastian Hoehna" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Func_time::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "Get the current system time in milliseconds.";
    
    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string Func_time::getHelpDetails(void) const
{
    // create a variable for the description of the function
    std::string details;
    return details;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Func_time::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    
    example += "time()\n";
    example += "\n";
    example += "# Wait a little bit\n";
    example += "sum = 0\n";
    example += "for (i in 1:10000) sum += i\n";
    example += "# Now print the time again\n";
    example += "time()\n";
    
    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> Func_time::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Func_time::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Func_time::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Get the time information";
    
    return title;
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
