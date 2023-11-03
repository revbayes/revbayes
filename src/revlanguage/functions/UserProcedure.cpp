#include <cstddef>
#include <sstream>
#include <list>
#include <vector>

#include "Signals.h"
#include "SyntaxElement.h"
#include "TypeSpec.h"
#include "UserFunctionDef.h"
#include "UserProcedure.h"
#include "Argument.h"
#include "Environment.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"

namespace RevLanguage { class ArgumentRules; }

using namespace RevLanguage;

/** Basic constructor */
UserProcedure::UserProcedure( UserFunctionDef* def ) :
    Function(),
    functionDef( def )
{
}


/** Clone the object */
UserProcedure* UserProcedure::clone(void) const
{
    return new UserProcedure(*this);
}


/**
 * Execute the procedure. We simply execute the code; a procedure is never part of
 * a deterministic variable, so we need not worry about dynamic, multiple execution
 * of the procedure from a deterministic variable.
 *
 * Arguments are added to our environment as aliases, so we can both modify them and
 * assign to them.
 */
RevPtr<RevVariable> UserProcedure::execute( void )
{
    // Clear signals
    Signals::getSignals().clearFlags();
    
    // Set initial return value
    RevPtr<RevVariable> retVar = NULL;
    
    // Create new evaluation frame with function base class execution environment as parent
    Environment* procedureFrame = new Environment( getEnvironment(), "UserProcedureEnvironment" );
    
    // Add the arguments to our environment as alias variables (modifiable and assignable)
    for (std::vector<Argument>::iterator it = args.begin(); it != args.end(); ++it)
    {
        procedureFrame->addAlias( it->getLabel(), it->getVariable() );
    }
    
    // Execute code
    const std::list<SyntaxElement*>& code = functionDef->getCode();
    for ( std::list<SyntaxElement*>::const_iterator i=code.begin(); i!=code.end(); i++ ) {
        
        SyntaxElement* theSyntaxElement = *i;
        retVar = theSyntaxElement->evaluateContent( *procedureFrame );
        
        if ( Signals::getSignals().isSet( Signals::RETURN ) )
        {
            Signals::getSignals().clearFlags();
            break;
        }
    }
    
    // Delete evaluation frame
    delete procedureFrame;
    
    // Return the return value
    return retVar;
}


/** Get Rev type (static) */
const std::string& UserProcedure::getClassType(void)
{
    static std::string rev_type = "UserProcedure";
    
	return rev_type;
}


/** Get Rev type spec (static) */
const TypeSpec& UserProcedure::getClassTypeSpec(void)
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), &Function::getClassTypeSpec() );
    
	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string UserProcedure::getFunctionName( void ) const
{
    // create a name variable that is NOT the same for all instance of this class
    std::string f_name = "UserProcedure";
    
    return f_name;
}


/** Get Rev type spec (from an instance) */
const TypeSpec& UserProcedure::getTypeSpec( void ) const
{
    return getClassTypeSpec();
}


/** Get argument rules */
const ArgumentRules& UserProcedure::getArgumentRules(void) const
{
    return functionDef->getArgumentRules();
}


/** Get return type */
const TypeSpec& UserProcedure::getReturnType(void) const
{
    return functionDef->getReturnType();
}

