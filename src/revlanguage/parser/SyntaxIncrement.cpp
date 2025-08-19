#include <cstddef>
#include <iostream>
#include <set>

#include "RbException.h"
#include "SyntaxIncrement.h"
#include "Environment.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "SyntaxElement.h"

using namespace RevLanguage;


/**
 * Construct from variable expression.
 *
 * \param[in]   v   The variable expression
 */
SyntaxIncrement::SyntaxIncrement( SyntaxElement* v, bool post ) : SyntaxElement(),
    postIncrement( post ),
    variable( v )
{
}


/**
 * Deep copy constructor. It really should not be possible for the
 * parser to return a NULL variable statement, so we do not have
 * to check for a NULL pointer.
 */
SyntaxIncrement::SyntaxIncrement( const SyntaxIncrement& x ) : SyntaxElement(x),
    postIncrement( x.postIncrement )
{
    variable   = x.variable->clone();
}


/**
 * Destructor deletes operands
 */
SyntaxIncrement::~SyntaxIncrement()
{
    delete variable;
}


/**
 * Assignment operator performing deep assignment.
 */
SyntaxIncrement& SyntaxIncrement::operator=( const SyntaxIncrement& x )
{
    if ( this != &x )
    {
        SyntaxElement::operator=( x );
        
        postIncrement = x.postIncrement;
        
        delete variable;
        variable = x.variable->clone();
    }
    
    return ( *this );
}


/**
 * The clone function is a convenience function to create proper copies of inherited objects.
 * That is, a.clone() will create a clone of type 'b' if the 'a' instance is of derived type 'b'.
 *
 * \return A new copy of myself
 */
SyntaxIncrement* SyntaxIncrement::clone () const
{
    return new SyntaxIncrement( *this );
}


/**
 * Evaluate the content of this syntax element. This will perform
 * an increment assignment operation.
 */
RevPtr<RevVariable> SyntaxIncrement::evaluateContent( const std::shared_ptr<Environment>& env, bool dynamic )
{
    
    RevPtr<RevVariable> retVar;
    
    // Get variable. We use standard evaluation because the variable is
    // implicitly on both sides (lhs and rhs) of this type of statement
    RevPtr<RevVariable> the_variable = variable->evaluateContent( env );
    if ( the_variable == NULL )
        throw RbException( "Invalid NULL variable returned by variable expression in increment assignment" );
    
    // Make sure that the variable is constant
    if ( !the_variable->getRevObject().isConstant() )
        throw RbException( "Invalid increment assignment to dynamic variable" );
    
    // Get a non-const reference to the lhs value object
    RevObject& lhs_value = the_variable->getRevObject();
    
    if ( postIncrement )
    {
        retVar = new RevVariable( lhs_value.clone() );
    }
    
    // Increment the lhs value. This will not change the control variable status.
    lhs_value.increment();
    
    if ( !postIncrement )
    {
        retVar = new RevVariable( lhs_value.clone() );
    }
        
    // No further assignment with this type of statement
    return retVar;
}


/**
 * Is the syntax element safe for use in a function (as
 * opposed to a procedure)? The assignment is safe
 * if the variable expression is safe, and the
 * increment is not to an external variable.
 */
bool SyntaxIncrement::isFunctionSafe( const Environment& env, std::set<std::string>& localVars ) const
{
    // Check variable
    if ( !variable->isFunctionSafe( env, localVars ) )
        return false;
    
    // Check whether assignment is to external variable (not function safe)
    // We are not allowed to add the variable to the local variables, it must
    // exist before the increment is attempted
    if ( variable->retrievesExternVar( env, localVars, false ) )
        return false;
    
    // All tests passed
    return true;
}


