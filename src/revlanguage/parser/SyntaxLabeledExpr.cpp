#include <cstddef>
#include <sstream>
#include <set>

#include "Environment.h"
#include "SyntaxLabeledExpr.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "SyntaxElement.h"

using namespace RevLanguage;

/** Constructor from label and expression */
SyntaxLabeledExpr::SyntaxLabeledExpr( const std::string& lbl, SyntaxElement* expr ) :
    SyntaxElement(),
    label( lbl ),
    expression( expr )
{
}


/** Deep copy constructor */
SyntaxLabeledExpr::SyntaxLabeledExpr( const SyntaxLabeledExpr& x ) :
SyntaxElement(x)
{
    label      = x.label;
    expression = x.expression->clone();
}


/** Destructor deletes members */
SyntaxLabeledExpr::~SyntaxLabeledExpr( void )
{
    delete expression;
}


/** Assignment operator */
SyntaxLabeledExpr& SyntaxLabeledExpr::operator=( const SyntaxLabeledExpr& x )
{
    if ( &x != this )
    {
        SyntaxElement::operator=(x);
        
        delete expression;

        label      = x.label;
        expression = x.expression->clone();
    }

    return (*this);
}


/** Clone syntax element */
SyntaxLabeledExpr* SyntaxLabeledExpr::clone( void ) const
{
    return new SyntaxLabeledExpr(*this);
}


/** Get semantic value (not applicable so return NULL) */
RevPtr<RevVariable> SyntaxLabeledExpr::evaluateContent( const std::shared_ptr<Environment>& env, bool dynamic )
{
    return NULL;
}


/** Is the expression constant?
 *  Only if the argument is constant.
 */
bool SyntaxLabeledExpr::isConstExpression( void ) const
{
    return expression->isConstExpression();
}


/**
 * Is this syntax element safe for inclusion in a function definition? Simply
 * delegate to the expression.
 */
bool SyntaxLabeledExpr::isFunctionSafe( const Environment& env, std::set<std::string>& localVars ) const
{
    return expression->isFunctionSafe( env, localVars );
}


/**
 * Does this syntax element retrieve an external variable? Simply
 * delegate to the expression. The inLHS flag should always be
 * false if we reach this element, and we pass it on as false.
 * This is because procedures are not allowed in LHS expressions,
 * and only procedures need to use this call.
 */
bool SyntaxLabeledExpr::retrievesExternVar( const Environment& env, std::set<std::string>& localVars, bool inLHS ) const
{
    return expression->retrievesExternVar( env, localVars, false );
}

