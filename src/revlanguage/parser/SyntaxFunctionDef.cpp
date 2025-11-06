#include <cstddef>
#include <sstream>
#include <list>
#include <set>
#include <string>

#include "Environment.h"
#include "RbException.h"
#include "SyntaxFunctionDef.h"
#include "UserFunction.h"
#include "UserFunctionDef.h"
#include "RlUserInterface.h"
#include "UserProcedure.h"
#include "Workspace.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "FunctionTable.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlUtils.h"
#include "SyntaxElement.h"
#include "SyntaxFormal.h"
#include "TypeSpec.h"

namespace RevLanguage { class Function; }

using namespace RevLanguage;

/** Construct global function call from function name and arguments */
SyntaxFunctionDef::SyntaxFunctionDef( const std::string&            type,
                                      const std::string&            name,
                                      std::list<SyntaxFormal*>*     formals,
                                      std::list<SyntaxElement*>*    stmts,
                                      bool                          isProcDef) :
    SyntaxElement(),
    code( stmts ),
    function_name( name ),
    formal_args( formals ),
    return_type( RlUtils::Void ),
    is_procedure_def( isProcDef )
{
    if ( type == "" )
        return_type = RevObject::getClassTypeSpec();
    else
        return_type  = TypeSpec( Workspace::userWorkspace().getClassTypeSpecOfType( type ) );
}


/** Deep copy constructor */
SyntaxFunctionDef::SyntaxFunctionDef( const SyntaxFunctionDef& x ) :
    SyntaxElement( x ),
    function_name( x.function_name ),
    return_type( x.return_type ),
    is_procedure_def( x.is_procedure_def )
{
    for ( std::list<SyntaxFormal*>::const_iterator it = x.formal_args->begin(); it != x.formal_args->end(); ++it )
        formal_args->push_back( (*it)->clone() );
    
    for ( std::list<SyntaxElement*>::const_iterator it = x.code->begin(); it != x.code->end(); ++it )
        code->push_back( (*it)->clone() );
}


/** Destructor deletes members */
SyntaxFunctionDef::~SyntaxFunctionDef( void )
{
    for ( std::list<SyntaxFormal*>::iterator it = formal_args->begin(); it != formal_args->end(); ++it )
        delete *it;
    delete formal_args;
    
    for ( std::list<SyntaxElement*>::iterator it = code->begin(); it != code->end(); ++it )
        delete *it;
    delete code;
}


/** Assignment operator */
SyntaxFunctionDef& SyntaxFunctionDef::operator=( const SyntaxFunctionDef& x )
{
    if ( &x != this )
    {
        SyntaxElement::operator=( x );

        function_name    = x.function_name;
        is_procedure_def  = x.is_procedure_def;
        return_type      = x.return_type;
        
        formal_args->clear();
        for ( std::list<SyntaxFormal*>::const_iterator it = x.formal_args->begin(); it != x.formal_args->end(); ++it )
            formal_args->push_back( (*it)->clone() );
        
        code->clear();
        for ( std::list<SyntaxElement*>::const_iterator it = x.code->begin(); it != x.code->end(); ++it )
            code->push_back( (*it)->clone() );
        
    }

    return ( *this );
}


/** Type-safe clone of syntax element */
SyntaxFunctionDef* SyntaxFunctionDef::clone( void ) const
{
    return new SyntaxFunctionDef( *this );
}



/**
 * Get semantic value: insert a user-defined function or procedure in the user workspace.
 *
 * @todo Deal with local variables hiding external variables. Ask if user wants to replace
 *       an existing function or procedure.
 */
RevPtr<RevVariable> SyntaxFunctionDef::evaluateContent( const std::shared_ptr<Environment>& env, bool dynamic )
{
    // Get argument rules from the formals
    ArgumentRules* argRules = new ArgumentRules();
    for ( std::list<SyntaxFormal*>::iterator it = formal_args->begin(); it !=formal_args->end(); ++it )
        argRules->push_back( (*it)->getArgumentRule()->clone() );

    // Check whether statements are function-safe if we are a function
    if ( !is_procedure_def )
    {
        // Load labeled formals as local variables
        std::set<std::string> localVars;
        for ( ArgumentRules::iterator it = argRules->begin(); it != argRules->end(); ++it )
        {
            if ( (*it)->getArgumentLabel() != "" )
                localVars.insert( (*it)->getArgumentLabel() );
        }

        // Now check that all statements are function-safe
        for( std::list<SyntaxElement*>::const_iterator it = code->begin(); it != code->end(); ++it )
        {
            if ( !(*it)->isFunctionSafe( *env, localVars ) )
            {
                std::ostringstream msg;
                msg << "The code of the function includes statements that modify or potentially modify" << std::endl;
                msg << "external variables. This is not allowed in a function. Either modify the code" << std::endl;
                msg << "of your function, or redefine it as a procedure." << std::endl;
                throw RbException( msg.str() );
            }
        }

        // Finally check whether last statement (if there is one) retrieves an external variable
        std::list<SyntaxElement*>::const_reverse_iterator rit = code->rbegin();
       if ( rit != code->rend() && (*rit)->retrievesExternVar( *env, localVars, false ) )
            throw RbException( "The code is not function-safe." );
    }

    // Create copy of the statements
    std::list<SyntaxElement*>* stmts = new std::list<SyntaxElement*>();

    for( std::list<SyntaxElement*>::const_iterator it = code->begin(); it != code->end(); ++it ) {
        stmts->push_back( (*it)->clone() );
    }
    
    // Create the function definition
    UserFunctionDef* functionDef = new UserFunctionDef( argRules, return_type, stmts, function_name );
    
    // Create the function or the procedure
    Function* the_function;
    if ( is_procedure_def )
        the_function = new UserProcedure( functionDef );
    else
        the_function = new UserFunction( functionDef );
    
    // Insert the function/procedure in the (user) workspace
    if ( env->getFunctionTable().existsFunctionInFrame( function_name, *argRules ) )
    {
        bool ok;
        if ( is_procedure_def )
            ok = UserInterface::userInterface().ask( "Replace existing procedure with same signature" );
        else
            ok = UserInterface::userInterface().ask( "Replace existing function with same signature" );

        if ( ok )
        {
            env->getFunctionTable().replaceFunction( function_name, the_function );
        }
        else
        {
            if ( is_procedure_def )
                RBOUT( "Registering of procedure canceled" );
            else
                RBOUT( "Registering of function canceled" );
        }
    }
    else
    {
        env->addFunction( the_function );
    }

    // No return value 
    return NULL;
}


