#include <cstddef>
#include <iostream>
#include <vector>
#include <set>
#include <string>

#include "Argument.h"
#include "Environment.h"
#include "RlFunction.h"
#include "SyntaxUnaryExpr.h"
#include "Workspace.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "SyntaxElement.h"

using namespace RevLanguage;


/** Static vector of string giving names of operator types */
std::string SyntaxUnaryExpr::opCode[] = { "uminus", "uplus", "unot" };


/** Construct from operator type and operands */
SyntaxUnaryExpr::SyntaxUnaryExpr( operatorT op, SyntaxElement* oper )  :
    SyntaxElement(),
    operand( oper ),
    operation(op)
{
}


/** Deep copy constructor */
SyntaxUnaryExpr::SyntaxUnaryExpr( const SyntaxUnaryExpr& x ) :
    SyntaxElement(x)
{
    operand     = x.operand->clone();
    operation   = x.operation;
}


/** Destructor deletes operand */
SyntaxUnaryExpr::~SyntaxUnaryExpr( void )
{
    delete operand;
}


/** Assignment operator */
SyntaxUnaryExpr& SyntaxUnaryExpr::operator=( const SyntaxUnaryExpr& x )
{
    if ( &x != this )
    {
        SyntaxElement::operator=( x );
        
        delete operand;
        
        operand     = x.operand->clone();
        operation   = x.operation;
    }

    return ( *this );
}



/** Type-safe clone of syntax element */
SyntaxUnaryExpr* SyntaxUnaryExpr::clone () const
{
    return new SyntaxUnaryExpr(*this);
}


/**
 * @brief Get semantic value (static version)
 *
 * We simply look up the function and calculate the value. Note that
 * we look up the function in the global workspace, which means that
 * we do not allow users to override internal functions corresponding
 * to unary expressions.
 *
 * Currently, we call the standard execute() method of the function,
 * which is typically set up to produce a deterministic variable.
 * This is more appropriate behavior for the dynamic version of this
 * function.
 *
 * It is not necessary to clone the function, we simply clear out
 * the argument after we are done.
 *
 * Even in the dynamic version of this function, it is not necessary to
 * clone the function. We simply use it to generate a deterministic
 * variable after setting the appropriate argument. After executing
 * the function, we make sure to clear out the argument of the function,
 * so that it is ready for the next round of execution.
 *
 * Note that we call the static evaluate function of the argument,
 * and make sure we return a constant value.
 *
 * Unlike the static version of this function, we call the dynamic evaluate
 * function of the argument. We also return the return value as is, without
 * making it a constant value first.
 *
 * @todo Support this evaluation context better
 */
RevPtr<RevVariable> SyntaxUnaryExpr::evaluateContent( const std::shared_ptr<Environment>& env, bool dynamic )
{
    // Package the argument
    std::vector<Argument> arg;
    arg.push_back( Argument( operand->evaluateContent( env, dynamic ), "" ) );
    
    // Find the function
    std::string funcName = "_" + opCode[ operation ];
    Function* func = Workspace::globalWorkspace().getFunction( funcName, arg, !dynamic ).clone();
    func->processArguments( arg );
    
    RevPtr<RevVariable> funcReturnValue = func->execute();
    
    // free the memory of our copy
    delete func;
    
    if ( dynamic == false || isConstExpression() == true )
    {
        // Return value after making sure it is constant
        if ( funcReturnValue != NULL )
        {
            funcReturnValue->getRevObject().makeConstantValue();
        }
        
    }
    
    return funcReturnValue;
}


/**
 * Is the expression constant?
 * Only if the argument (operand) is constant.
 */
bool SyntaxUnaryExpr::isConstExpression( void ) const
{
    return operand->isConstExpression();
}


/**
 * Is the syntax element safe for use in a function (as
 * opposed to a procedure)? The unary expression is safe
 * if its single operand is safe.
 */
bool SyntaxUnaryExpr::isFunctionSafe( const Environment& env, std::set<std::string>& localVars ) const
{
    // Check the operand
    if ( operand->isFunctionSafe( env, localVars ) )
        return true;
    else
        return false;
}

