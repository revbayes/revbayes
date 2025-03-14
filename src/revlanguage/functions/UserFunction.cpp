#include <cstddef>
#include <sstream>
#include <list>
#include <string>
#include <vector>
#include <memory>

#include "RbException.h"
#include "Signals.h"
#include "TypeSpec.h"
#include "UserFunction.h"
#include "UserFunctionDef.h"
#include "Workspace.h"
#include "Argument.h"
#include "Environment.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "SyntaxElement.h"

namespace RevBayesCore { class DagNode; }
namespace RevLanguage { class ArgumentRules; }

using namespace RevLanguage;

/** Basic constructor */
UserFunction::UserFunction( UserFunctionDef* def ) : Function(),
    function_def( def )
{
}


/** Clone the object */
UserFunction* UserFunction::clone(void) const
{
    return new UserFunction(*this);
}


/** Execute function. Here we create a deterministic node if applicable, otherwise we just execute the code */
RevPtr<RevVariable> UserFunction::execute( void )
{
    
    // If the return type object has a DAG node inside it, we return an appropriate model/container/factor object
    // with a deterministic node inside it. Otherwise we return a "flat" RevObject without a dag node inside it.

    RevObject* retVal = Workspace::userWorkspace().makeNewDefaultObject( getReturnType().getType() );

    if ( retVal->isModelObject() )
    {
        retVal->makeUserFunctionValue( this->clone() );

        return new RevVariable( retVal );
    }
    else
    {
        // "Flat" call: Simply execute and return the variable
        delete retVal;  // We don't need this

        return executeCode();
    }
}


/** In this function we execute the Rev code for the function (compiled syntax tree) */
RevPtr<RevVariable> UserFunction::executeCode( void )
{
    // Create new evaluation frame with function base class execution environment as parent
    auto function_frame = std::make_shared<Environment>( getEnvironment(), "UserFunctionEnvironment" );
    
    // Add the arguments to our environment
    for ( std::vector<Argument>::iterator it = args.begin(); it != args.end(); ++it )
    {
        // Note: We can add also temporary variable arguments as references because we
        // currently store them as arguments of the Rev function in UserFunctionArgs
        // as long as the UserFunctionCall exists.
        function_frame->addReference( it->getLabel(), it->getVariable() );
    }

    // Clear signals
    Signals::getSignals().clearFlags();
    
    // Set initial return value
    RevPtr<RevVariable> ret_var = NULL;
    
    // Execute code
    const std::list<SyntaxElement*>& code = function_def->getCode();
    for ( std::list<SyntaxElement*>::const_iterator it = code.begin(); it != code.end(); ++it )
    {
        SyntaxElement* the_syntax_element = *it;
        ret_var = the_syntax_element->evaluateContent( function_frame );
        
        if ( Signals::getSignals().isSet( Signals::RETURN ) )
        {
            Signals::getSignals().clearFlags();
            break;
        }
    }

    // non-void return type?
    if ( getReturnType() != RevObject::getClassTypeSpec() )
    {
        // void return value?
        if (ret_var == NULL)
        {
            throw(RbException("No return value in function '"+this->getFunctionName()+"' returning non-void type "+getReturnType().getType()));
        }
        else if ( ret_var->getRevObject().isType( getReturnType() ) == true )
        {
            if ( ret_var->getRevObject().getTypeSpec() != getReturnType() && ret_var->getRevObject().getTypeSpec().isDerivedOf( getReturnType() ) == false )
            {
                // compatible but differing return value
                if ( ret_var->getRevObject().isConvertibleTo(getReturnType(),true) >= 0 )
                {
                    //convert the return value
                    ret_var = new RevVariable( ret_var->getRevObject().convertTo(getReturnType()) );
                }
                // incompatible return value?
                else
                {
                    throw(RbException("Returning "+ret_var->getRevObject().getTypeSpec().getType()+" in function '"+this->getFunctionName()+"' with incompatible return type "+getReturnType().getType()));
                }
            }
        }
        else
        {
            throw(RbException("Returning "+ret_var->getRevObject().getTypeSpec().getType()+" in function '"+this->getFunctionName()+"' with incompatible return type "+getReturnType().getType()));
        }
    }

    // Return the return value
    return ret_var;
}


/** Get Rev type (static) */
const std::string& UserFunction::getClassType(void)
{
    static std::string rev_type = "UserFunction";
    
	return rev_type; 
}


/** Get Rev type spec (static) */
const TypeSpec& UserFunction::getClassTypeSpec(void)
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), &Function::getClassTypeSpec() );
    
	return rev_type_spec; 
}


/** Get the parameters from the argument vector */
std::vector<const RevBayesCore::DagNode*> UserFunction::getParameters(void) const
{
    std::vector<const RevBayesCore::DagNode*> params;

    for (std::vector<Argument>::const_iterator it = args.begin(); it != args.end(); ++it )
    {
        
        if ( (*it).getVariable()->getRevObject().isModelObject() )
        {
            params.push_back( (*it).getVariable()->getRevObject().getDagNode() );
        }
        
    }

    return params;
}


/**
 * Get the primary Rev name for this function.
 */
std::string UserFunction::getFunctionName( void ) const
{
    // create a name variable that is NOT the same for all instance of this class
    std::string f_name = function_def->getName();
    
    return f_name;
}


/** Get Rev type spec (from an instance) */
const TypeSpec& UserFunction::getTypeSpec( void ) const
{
    return getClassTypeSpec();
}


/** Get argument rules */
const ArgumentRules& UserFunction::getArgumentRules(void) const
{
    return function_def->getArgumentRules();
}


/** Get return type */
const TypeSpec& UserFunction::getReturnType(void) const
{
    return function_def->getReturnType();
}

