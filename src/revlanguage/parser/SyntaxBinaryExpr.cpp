#include <cstddef>
#include <iostream>
#include <vector>
#include <set>

#include "Argument.h"
#include "RlFunction.h"
#include "SyntaxBinaryExpr.h"
#include "Workspace.h"
#include "Environment.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "SyntaxElement.h"

using namespace RevLanguage;

/** Static vector of strings giving names of operator types */
std::string SyntaxBinaryExpr::opCode[] = { "range", "add", "sub", "mul", "div", "mod", "power", "lt", "le",
                                           "eq", "ne", "ge", "gt", "and", "or", "and", "or"};


/** Construct from operator type and operands */
SyntaxBinaryExpr::SyntaxBinaryExpr( operatorT op, SyntaxElement* lhs, SyntaxElement* rhs ) :
    SyntaxElement(),
    left_operand( lhs ),
    right_operand( rhs ),
    operation( op )
{
}


/** Deep copy constructor */
SyntaxBinaryExpr::SyntaxBinaryExpr( const SyntaxBinaryExpr& x ) :
    SyntaxElement(x)
{
    left_operand  = x.left_operand->clone();
    right_operand = x.right_operand->clone();
    operation    = x.operation;
}


/** Destructor deletes operands */
SyntaxBinaryExpr::~SyntaxBinaryExpr()
{
    delete left_operand;
    delete right_operand;
}


/** Assignment operator */
SyntaxBinaryExpr& SyntaxBinaryExpr::operator=( const SyntaxBinaryExpr& x )
{

    if ( &x != this )
    {
        SyntaxElement::operator=( x );
        
        delete left_operand;
        delete right_operand;

        left_operand  = x.left_operand->clone();
        right_operand = x.right_operand->clone();
        operation    = x.operation;
    }

    return *this;
}


/** Type-safe clone of syntax element */
SyntaxElement* SyntaxBinaryExpr::clone () const
{
    return (new SyntaxBinaryExpr(*this));
}


/**
 * @brief Get semantic value
 *
 * We simply look up the function and calculate the value. Note that
 * we look up the function in the global workspace, which means that
 * we do not allow users to override internal functions corresponding
 * to binary expressions.
 *
 * Currently, we call the standard execute() method of the function,
 * which is typically set up to produce a deterministic variable.
 * This is more appropriate behavior for the dynamic version of this
 * function.
 *
 * It is not necessary to clone the function, we simply clear out
 * the arguments after we are done.
 *
 * Even in the dynamic version of this function, it is not necessary to
 * clone the function. We simply use it to generate a deterministic
 * variable after setting the appropriate arguments. After executing
 * the function, we make sure to clear out the arguments of the function,
 * so that it is ready for the next round of execution.
 *
 * If dynamic == false, then we call the static evaluate function of the arguments,
 * and make sure we return a constant value.
 *
 * If dynamic == true, then we call the dynamic evaluate
 * function of the arguments. We also return the return value as is, without
 * making it a constant value first.
 *
 * @todo Support this evaluation context better
 */
RevPtr<RevVariable> SyntaxBinaryExpr::evaluateContent( const std::shared_ptr<Environment>& env, bool dynamic )
{
    
    // Package the arguments
    std::vector<Argument> args;
    
    RevPtr<RevVariable> left = left_operand->evaluateContent( env, dynamic );
    args.push_back( Argument( left, "" ) );
    
    RevPtr<RevVariable> right = right_operand->evaluateContent( env, dynamic );
    args.push_back( Argument( right, "" ) );
    
    // Get function and create deterministic DAG node
    std::string func_name = "";
    
    if ( operation != Range && operation != Exp)
    {
        func_name += "_";
    }
    
    func_name += opCode[ operation ];
    Function* the_function = Workspace::globalWorkspace().getFunction( func_name, args, false ).clone();
    the_function->processArguments( args );
    
    RevPtr<RevVariable> the_return_value = the_function->execute();
    
    // Free the memory of our copy
    delete the_function;
    
    if ( isConstExpression() == true )
    {
        // Return the return value of the function after making it constant
        if ( the_return_value != NULL )
        {
            the_return_value->getRevObject().makeConstantValue();
        }
        
    }
    return the_return_value;
}


/** Is the expression constant?
 *  Only if the arguments are constant.
 */
bool SyntaxBinaryExpr::isConstExpression( void ) const
{
    return left_operand->isConstExpression() && right_operand->isConstExpression();
}


/**
 * Is the syntax element safe for use in a function (as
 * opposed to a procedure)? The binary expression is safe
 * if its operands are safe.
 */
bool SyntaxBinaryExpr::isFunctionSafe( const Environment& env, std::set<std::string>& localVars ) const
{
    // Check operands
    if ( left_operand->isFunctionSafe( env, localVars ) && right_operand->isFunctionSafe( env, localVars ) )
        return true;

    // At least one operand not safe
    return false;
}

