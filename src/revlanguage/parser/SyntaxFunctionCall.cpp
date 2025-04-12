#include <cstddef>
#include <sstream>
#include <list>
#include <set>
#include <vector>

#include "Argument.h"
#include "Environment.h"
#include "RlMemberMethod.h"
#include "RbException.h"
#include "RlString.h"
#include "SyntaxFunctionCall.h"
#include "Workspace.h"
#include "MethodTable.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "SyntaxElement.h"
#include "SyntaxPipePlaceholder.h"
#include "SyntaxLabeledExpr.h"
#include "SyntaxVariable.h"

using namespace RevLanguage;

/** Construct global function call from function name and arguments */
SyntaxFunctionCall::SyntaxFunctionCall( const std::string &n, std::list<SyntaxLabeledExpr*>* args) :
    SyntaxElement(),
    arguments( args ),
    function_name( n ),
    base_variable( NULL )
{
}


/** Construct member function call from base variable, function name and arguments */
SyntaxFunctionCall::SyntaxFunctionCall( SyntaxVariable* var, const std::string& n, std::list<SyntaxLabeledExpr*>* args ) :
    SyntaxElement(),
    arguments( args ),
    function_name( n ),
    base_variable( var )
{
}


/** Deep copy constructor */
SyntaxFunctionCall::SyntaxFunctionCall(const SyntaxFunctionCall& x) :
    SyntaxElement( x ),
    arguments( NULL ),
    function_name( x.function_name ),
    base_variable( NULL )
{
    if (x.base_variable != NULL)
        base_variable = x.base_variable->clone();
    
    arguments = new std::list<SyntaxLabeledExpr*>();
    for ( std::list<SyntaxLabeledExpr*>::iterator it = x.arguments->begin(); it != x.arguments->end(); ++it )
        arguments->push_back( (*it)->clone() );
}


/** Destructor deletes members */
SyntaxFunctionCall::~SyntaxFunctionCall()
{
    delete base_variable;
    
    for ( std::list<SyntaxLabeledExpr*>::iterator it = arguments->begin(); it != arguments->end(); ++it )
        delete *it;
    delete arguments;
}


/** Assignment operator */
SyntaxFunctionCall& SyntaxFunctionCall::operator=( const SyntaxFunctionCall& x )
{
    if ( &x != this )
    {
        SyntaxElement::operator=( x );

        delete base_variable;
        
        for ( std::list<SyntaxLabeledExpr*>::iterator it = arguments->begin(); it != arguments->end(); ++it )
        {
            delete *it;
        }
        
        function_name = x.function_name;

        if (x.base_variable != NULL)
        {
            base_variable = x.base_variable->clone();
        }
        
        arguments->clear();
        for ( std::list<SyntaxLabeledExpr*>::iterator it = x.arguments->begin(); it != x.arguments->end(); ++it )
        {
            arguments->push_back( (*it)->clone() );
        }
        
    }

    return ( *this );
}


/** Type-safe clone of syntax element */
SyntaxFunctionCall* SyntaxFunctionCall::clone( void ) const
{
    return new SyntaxFunctionCall( *this );
}



/**
 * Evaluate semantic content.
 *
 * If dynamic == false, then we know that we only want to evaluate the function once,
 * and get a constant value back. For now, it does the same thing as
 * if dynamic == true, with two exceptions. First, it uses the evaluate function
 * with dynamic == false to evaluate its arguments. Second, it makes the return
 * value a constant value.
 *
 * If dynamic == true, then we know that we want
 * to evaluate the function repeatedly, as part of a deterministic variable.
 * Therefore, we need to make a deterministic variable with the function inside
 * it. Currently, this is the default implementation for the execute function
 * in RevLanguage::Function, but we may want to differentiate between the two
 * contexts later.
 *
 * @todo Support this function call context better
 */
RevPtr<RevVariable> SyntaxFunctionCall::evaluateContent( const std::shared_ptr<Environment>& env, bool dynamic )
{
    
    // Package arguments
    std::vector<Argument> args;
    for ( std::list<SyntaxLabeledExpr*>::const_iterator it = arguments->begin(); it != arguments->end(); ++it )
    {
        
        const RlString& theLabel = (*it)->getLabel();
        RevPtr<RevVariable> theVar = (*it)->getExpression().evaluateContent(env,dynamic);
        
        Argument theArg = Argument( theVar, theLabel.getValue() );
        args.push_back( theArg );
    }
    
    Function* func = NULL;
    if ( base_variable == NULL )
    {
        // We are trying to find a function in the current environment
        
        // First we see if the function name corresponds to a user-defined variable
        // We can do this first because user-defined variables are not allowed to mask function names
        // Skip if we're not in UserWorkspace, because functions can only be user-defined in UserWorkspace
        bool found = false;
        if ( env->existsVariable( function_name ) && env == Workspace::userWorkspacePtr() )
        {
            RevObject &the_object = env->getRevObject( function_name );
            
            if ( the_object.isType( Function::getClassTypeSpec() ) )
            {
                func = static_cast<Function*>( the_object.clone() );
                std::vector<bool> arg_mapped(args.size(), false);
                found = func->checkArguments(args, NULL, arg_mapped);
            }
        }
        
        // If we cannot find the function name as a variable, it must be in the function table
        // This call will throw a relevant message if the function is not found
        if ( found == false )
        {
            func = env->getFunction(function_name, args, !dynamic).clone();
        }
        
        // Allow the function to process the arguments
        func->processArguments( args );
        
        // Set the execution environment of the function
        func->setExecutionEnviroment( env );
    }
    else
    {
        // We are trying to find a member function
        
        // First we get the base variable
        RevPtr<RevVariable> the_var = base_variable->evaluateContent( env, dynamic );
        
        // Now we get a reference to the member object inside
        RevObject &the_member_object = the_var->getRevObject();
        
        const MethodTable& mt = the_member_object.getMethods();
        
        const Function* the_const_function = mt.findFunction( function_name, args, !dynamic );

        Function* the_function;
        if ( the_const_function != NULL )
        {
            the_function = the_const_function->clone();
        }
        else
        {
            throw RbException()<<"Variable of type '"<<the_member_object.getType()<<"' has no method called '"<<function_name<<"'.  You can use '.methods()' to find available methods.";
        }
        
        the_function->processArguments(args );
        
        MemberMethod* the_member_method = dynamic_cast<MemberMethod*>( the_function );
        if ( the_member_method != NULL )
        {
            the_member_method->setMemberObject( the_var );
            func = the_function;
        }
        else
        {
            delete the_function;
            throw RbException("Error while trying to access member function/procedure.");
        }
        
    }
    
    // Evaluate the function
    RevPtr<RevVariable> func_return_value = func->execute();
    
    // free the memory of our copy
    delete func;
    
    if ( isConstExpression() == true )
    {
        // Return the value, which is typically a deterministic variable with the function
        // inside it, although many functions return constant values or NULL (void).
        // To make sure the value is a constant and not a deterministic variable in this
        // context, we convert the return value here to a constant value before returning it.
        if ( func_return_value != NULL )
        {
            func_return_value->getRevObject().makeConstantValue();
        }
    }
    
    return func_return_value;
}

std::pair<int,int> SyntaxFunctionCall::pipeAddArgPlaceholder(SyntaxElement* piped_arg)
{
    assert(piped_arg);
    int n_placeholders = 0;

    // 1. Check arguments of final fxncall for _
    for(auto& argument: *arguments)
    {
        // If the argument is of the form `_` or `label = _`, the replace the `_` with piped_arg.
        if (dynamic_cast<SyntaxPipePlaceholder*>(&argument->getExpression()))
        {
            n_placeholders++;
            std::string label = argument->getLabel();
            delete argument;
            argument = new SyntaxLabeledExpr(label, piped_arg);
        }
    }

    int n_fxncalls = 1;
    if (base_variable)
    {
        // 2. Check if base_variable is _
        if (dynamic_cast<SyntaxPipePlaceholder*>(base_variable))
        {
            n_placeholders++;
            delete base_variable;
            base_variable = piped_arg;
        }
        // 3. Check if base_variable is another function call
        else if (auto sub_fxncall = dynamic_cast<SyntaxFunctionCall*>(base_variable))
        {
            auto tmp = sub_fxncall->pipeAddArgPlaceholder(piped_arg);
            n_placeholders += tmp.first;
            n_fxncalls += tmp.second;
        }
    }

    return {n_placeholders, n_fxncalls};
}

void SyntaxFunctionCall::pipeAddArg(SyntaxElement* piped_arg)
{
    auto tmp = pipeAddArgPlaceholder(piped_arg);
    int n_placeholders   = tmp.first;
    int n_fxncalls = tmp.second;

    if (n_placeholders > 1)
        throw RbException()<<"Pipe placeholder may only occur once.";

    if (n_placeholders == 0)
    {
        if (n_fxncalls > 1)
            throw RbException()<<"Piping into expression with multiple function calls requires a pipe placeholder ('_') !";

        arguments->push_front(new SyntaxLabeledExpr ("" , piped_arg));
    }
}

/**
 * Is the expression constant?
 * Only if all arguments are constant.
 */
bool SyntaxFunctionCall::isConstExpression(void) const
{
    
    // We need to iterate over all arguments
    for ( std::list<SyntaxLabeledExpr*>::const_iterator it = arguments->begin(); it != arguments->end(); ++it )
    {
        // We return false if this argument is not constant
        SyntaxLabeledExpr* expr = *it;
        if ( expr->isConstExpression() == false )
        {
            return false;
        }
    }
    
    if ( base_variable != NULL )
    {
        if ( base_variable->isConstExpression() == false )
        {
            return false;
        }
    }
    
    // All arguments are constant
    return true;
}


/**
 * Is the syntax element safe for use in a function
 * (as opposed to a procedure)? The function call is safe
 * if it is a call to a function, and the argument expressions
 * are all function-safe. If it is a call to a procedure, it
 * is safe if all argument expressions are function-safe, and
 * none of them retrieves an external variable.
 */
bool SyntaxFunctionCall::isFunctionSafe( const Environment& env, std::set<std::string>& localVars ) const
{
    // Protect from self-checking if recursive. If that case, the function
    // does not exist yet and we tentatively assume it is safe
    if ( !env.existsFunction( function_name ) )
        return true;
    
    if ( env.isProcedure( function_name ) )
    {
        // Check base variable
        if ( base_variable != NULL && ( !base_variable->isFunctionSafe( env, localVars ) || base_variable->SyntaxElement::retrievesExternVar( env, localVars, false ) ) )
            return false;
        
        // Iterate over all arguments
        for ( std::list<SyntaxLabeledExpr*>::const_iterator it = arguments->begin(); it != arguments->end(); ++it )
        {
            // Return false if argument expression is not function-safe or retrieves an external variable
            SyntaxLabeledExpr* expr = *it;
            if ( !expr->isFunctionSafe( env, localVars ) || expr->retrievesExternVar( env, localVars, false ) )
                return false;
        }
        
    }
    else
    {
        // Check base variable
        if ( base_variable != NULL && !base_variable->isFunctionSafe( env, localVars ) )
            return false;
        
        // Iterate over all arguments
        for ( std::list<SyntaxLabeledExpr*>::const_iterator it = arguments->begin(); it != arguments->end(); ++it )
        {
            // Return false if argument expression is not function-safe
            SyntaxLabeledExpr* expr = *it;
            if ( !expr->isFunctionSafe( env, localVars ) )
                return false;
        }
        
    }
    
    // All arguments and the base variable are OK
    return true;
}

