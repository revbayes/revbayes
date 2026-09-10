#include <iostream>
#include <set>

#include "SyntaxReferenceAssignment.h"
#include "Environment.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "SyntaxElement.h"

using namespace RevLanguage;

/** Basic constructor from lef-hand side and right-hand side expressions */
SyntaxReferenceAssignment::SyntaxReferenceAssignment( SyntaxElement* lhsExpr, SyntaxElement* rhsExpr ) :
    SyntaxElement(),
    lhsExpression( lhsExpr ),
    rhsExpression( rhsExpr )
{
}


/** Deep copy constructor */
SyntaxReferenceAssignment::SyntaxReferenceAssignment( const SyntaxReferenceAssignment& x ) :
    SyntaxElement( x ),
    lhsExpression( x.lhsExpression->clone() ),
    rhsExpression( x.rhsExpression->clone() )
{
}


/** Destructor deletes operands */
SyntaxReferenceAssignment::~SyntaxReferenceAssignment( void )
{
    delete lhsExpression;
    delete rhsExpression;
}


/** Assignment operator */
SyntaxReferenceAssignment& SyntaxReferenceAssignment::operator=( const SyntaxReferenceAssignment& x )
{
    if ( this != &x ) {
        
        delete lhsExpression;
        delete rhsExpression;
        
        lhsExpression = x.lhsExpression->clone();
        rhsExpression = x.rhsExpression->clone();
    }
    
    return ( *this );
}


/** Type-safe clone of the syntax element */
SyntaxReferenceAssignment* SyntaxReferenceAssignment::clone () const
{
    return new SyntaxReferenceAssignment( *this );
}


/**
 * Get semantic value. When evaluating the semantic value of a reference assignment,
 * we first evaluate the rhs expression as if it were a constant expression. Among
 * other things, this makes it possible for us to create references to control
 * variables, which would be impossible if the rhs expression was evaluated as a
 * dynamic expression. We are interested in creating a reference to the variable
 * that results from evaluation of the rhs expression now.
 *
 * Note that the return variable is variable returned by the rhs expression.
 * We need not clone it.
 */
RevPtr<RevVariable> SyntaxReferenceAssignment::evaluateContent( const std::shared_ptr<Environment>& env, bool dynamic )
{
    
    // Declare variable storing the return value of the assignment expression
    RevPtr<RevVariable> the_variable;
    
    // Get the rhs expression wrapped and executed into a variable.
    the_variable = rhsExpression->evaluateContent( env );
    
    // Get variable slot from lhs
    RevPtr<RevVariable> theSlot;
    theSlot = lhsExpression->evaluateLHSContent( env, the_variable->getRevObject().getType() );
    
    // Make the slot a reference to the rhs expression variable.
    theSlot->makeReference( the_variable );
    
    // Return variable
    return the_variable;
}


/** This is an assignment, return true. */
bool SyntaxReferenceAssignment::isAssignment( void ) const
{
    return true;
}


/**
 * Is the syntax element safe for use in a function (as
 * opposed to a procedure)? The assignment is safe
 * if its lhs and rhs expressions are safe, and the
 * assignment is not to an external variable.
 */
bool SyntaxReferenceAssignment::isFunctionSafe( const Environment& env, std::set<std::string>& localVars ) const
{
    // Check lhs and rhs expressions
    if ( !lhsExpression->isFunctionSafe( env, localVars ) || !rhsExpression->isFunctionSafe( env, localVars ) )
        return false;
    
    // Check whether assignment is to external variable (not function-safe)
    if ( lhsExpression->retrievesExternVar( env, localVars, true ) )
        return false;
    
    // All tests passed
    return true;
}


