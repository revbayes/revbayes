#include <cstddef>
#include <iostream>
#include <set>

#include "RbException.h"
#include "SyntaxDecrement.h"
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
SyntaxDecrement::SyntaxDecrement( SyntaxElement* v, bool post ) : SyntaxElement(),
    postDecrement( post ),
    variable( v )
{
}


/** 
 * Deep copy constructor. It really should not be possible for the
 * parser to return a NULL variable statement, so we do not have
 * to check for a NULL pointer.
 */
SyntaxDecrement::SyntaxDecrement( const SyntaxDecrement& x ) : SyntaxElement(x),
    postDecrement( x.postDecrement )
{
    variable   = x.variable->clone();
}


/** 
 * Destructor deletes operands 
 */
SyntaxDecrement::~SyntaxDecrement() 
{
    delete variable;
}


/** 
 * Assignment operator performing deep assignment.
 */
SyntaxDecrement& SyntaxDecrement::operator=( const SyntaxDecrement& x ) 
{
    if ( this != &x )
    {
        SyntaxElement::operator=( x );

        postDecrement = x.postDecrement;
        
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
SyntaxDecrement* SyntaxDecrement::clone () const
{
    return new SyntaxDecrement( *this );
}


/** 
 * Evaluate the content of this syntax element. This will perform
 * a decrement assignment operation.
 */
RevPtr<RevVariable> SyntaxDecrement::evaluateContent( const std::shared_ptr<Environment>& env, bool dynamic )
{
    
    RevPtr<RevVariable> retVar;
    
    // Get variable. We use standard evaluation because the variable is
    // implicitly on both sides (lhs and rhs) of this type of statement
    RevPtr<RevVariable> the_variable = variable->evaluateContent( env );
    if ( the_variable == NULL )
        throw RbException( "Invalid NULL variable returned by variable expression in decrement assignment" );
    
    // Make sure that the variable is constant
    if ( !the_variable->getRevObject().isConstant() )
        throw RbException( "Invalid decrement assignment to dynamic variable" );
    
    // Get a non-const reference to the lhs value object
    RevObject& lhs_value = the_variable->getRevObject();
    
    if ( postDecrement )
    {
        retVar = new RevVariable( lhs_value.clone() );
    }

    // Decrement the lhs value. This will not change the control variable status.
    lhs_value.decrement();
    
    if ( !postDecrement )
    {
        retVar = new RevVariable( lhs_value.clone() );
    }
    
    // Return the variable for further assignment
    return retVar;
}



/**
 * Is the syntax element safe for use in a function (as
 * opposed to a procedure)? The assignment is safe
 * if the variable expression is safe, and the
 * decrement is not to an external variable.
 */
bool SyntaxDecrement::isFunctionSafe( const Environment& env, std::set<std::string>& localVars ) const
{
    // Check variable
    if ( !variable->isFunctionSafe( env, localVars ) )
        return false;
    
    // Check whether assignment is to external variable (not function safe)
    // We are not allowed to add the variable to the local variables, it must
    // exist before the decrement is attempted
    if ( variable->retrievesExternVar( env, localVars, false ) )
        return false;
    
    // All tests passed
    return true;
}

