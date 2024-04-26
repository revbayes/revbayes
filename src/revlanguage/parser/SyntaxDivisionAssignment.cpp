#include <cstddef>

#include "RbException.h"
#include "SyntaxDivisionAssignment.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "SyntaxAssignment.h"

namespace RevLanguage { class SyntaxElement; }

using namespace RevLanguage;


/**
 * Construct from lhs and rhs expressions.
 * Here only the object is instantiated and the variables set.
 *
 * \param[in]   lhsExpr     The left-hand side expression.
 * \param[in]   rhsExpr     The right-hand side expression.
 */
SyntaxDivisionAssignment::SyntaxDivisionAssignment( SyntaxElement* lhsExpr, SyntaxElement* rhsExpr ) : SyntaxAssignment(lhsExpr,rhsExpr)
{
    
}


/**
 * Destructor deletes operands
 */
SyntaxDivisionAssignment::~SyntaxDivisionAssignment()
{

}


/**
 * The clone function is a convenience function to create proper copies of inherited objects.
 * That is, a.clone() will create a clone of type 'b' if the 'a' instance is of derived type 'b'.
 *
 * \return A new copy of myself
 */
SyntaxDivisionAssignment* SyntaxDivisionAssignment::clone () const
{
    return new SyntaxDivisionAssignment( *this );
}


/**
 * Evaluate the content of this syntax element. This will perform a
 * division assignment operation.
 */
void SyntaxDivisionAssignment::assign(RevPtr<RevVariable> &lhs, RevPtr<RevVariable> &rhs)
{
    
    // Get variable from lhs. We use standard evaluation because the variable is
    // implicitly on both sides (lhs and rhs) of this type of statement
    if ( lhs == NULL )
    {
        throw RbException( "Invalid NULL variable returned by lhs expression in division assignment" );
    }
    
    // Make sure that the variable is constant
    if ( !lhs->getRevObject().isConstant() )
    {
        throw RbException( "Invalid division assignment to dynamic variable" );
    }
    
    // Record whether it is a control variable
    bool isWorkspaceVar = lhs->isWorkspaceVariable();
    
    // Get a reference to the lhs value object
    const RevObject& lhs_value = lhs->getRevObject();
    
    // Evaluate the rhs expression
    if ( rhs == NULL )
    {
        throw RbException( "Invalid NULL variable returned by rhs expression in division assignment" );
    }
    
    // Get a reference to the rhs value object
    const RevObject& rhs_value = rhs->getRevObject();
    
    // Generate result of the multiplication
    RevObject *newValue = lhs_value.divide( rhs_value );
    
    // Fill the slot with the new RevVariable
    lhs->replaceRevObject( newValue );
    
    // Reset it as control variable, if it was a control variable before the assignment.
    // When we fill the slot, the control variable property is reset to false by default.
    if ( isWorkspaceVar )
    {
        lhs->setWorkspaceVariableState( true );
    }
    
}


