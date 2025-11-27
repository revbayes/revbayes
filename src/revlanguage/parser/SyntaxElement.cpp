#include <iosfwd>
#include <set>

#include "SyntaxElement.h"
#include "RevPtr.h"
#include "RevVariable.h"

namespace RevLanguage { class Environment; }

using namespace RevLanguage;


/**
 * This function evaluates the semantic value of the element if it occurs on the left-hand
 * side of an assignment. It allows variable syntax elements to do the right thing when they
 * appear on the left-hand side of an assignment. The default behavior is to evaluate the
 * semantic value in the standard way for constant expressions, which is appropriate for all
 * other types of syntax elements when they appear on the left-hand sign of an assignment.
 */
RevPtr<RevVariable> SyntaxElement::evaluateLHSContent( const std::shared_ptr<Environment>& env, const std::string& varType )
{
    return evaluateContent( env );
}


/**
 * Is the syntax element an assignment element?
 * By default the answer is no.
 */
bool SyntaxElement::isAssignment( void ) const
{
    return false;
}


/**
 * Is the syntax element a constant expression?
 * By default a syntax element is not constant.
 */
bool SyntaxElement::isConstExpression( void ) const
{
    return false;
}


/**
 * Is the syntax element safe for use in a function
 * (as opposed to a procedure)? Most elements are safe,
 * so we return true by default.
 */
bool SyntaxElement::isFunctionSafe( const Environment& env, std::set<std::string>& localVars ) const
{
    return true;
}


/**
 * Does the syntax element retireve an external variable?
 * This is a relevant question in the evaluation of whether
 * the syntax tree is safe for use in a function (as opposed
 * to a procedure). Most elements do not retrieve variables,
 * so we return false by default.
 */
bool SyntaxElement::retrievesExternVar( const Environment& env, std::set<std::string>& localVars, bool inLHS ) const
{
    return false;
}

