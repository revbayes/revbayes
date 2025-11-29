#include <cstddef>
#include <cassert>
#include <sstream>
#include <set>

#include "RlContainer.h"
#include "Environment.h"
#include "RbException.h"
#include "SyntaxForLoop.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "SyntaxElement.h"

using namespace RevLanguage;


/** Standard constructor */
SyntaxForLoop::SyntaxForLoop( const std::string& identifier, SyntaxElement* inExpr ) : SyntaxElement(),
    varName( identifier ),
    inExpression( inExpr ),
    stateSpace( NULL ),
    nextIndex( 0 )
{
    if ( inExpression == NULL )
    {
        throw RbException("The 'in' expression of for loop is empty");
    }
    
}


/** Deep copy constructor */
SyntaxForLoop::SyntaxForLoop( const SyntaxForLoop& x ) : SyntaxElement(x)
{
    varName                     = x.varName;
    inExpression                = x.inExpression->clone();
    stateSpace                  = NULL;
    nextIndex                   = 0;
}


/** Destructor deletes members */
SyntaxForLoop::~SyntaxForLoop()
{
    delete inExpression;
    delete stateSpace;
}


/** Assignment operator */
SyntaxForLoop& SyntaxForLoop::operator=( const SyntaxForLoop& x )
{
    if ( &x != this )
    {
        SyntaxElement::operator=(x);

        delete inExpression;
        delete stateSpace;

        varName                     = x.varName;
        inExpression                = x.inExpression->clone();
        stateSpace                  = NULL;
        nextIndex                   = 0;
    }

    return *this;
}


/** Clone syntax element */
SyntaxElement* SyntaxForLoop::clone( void ) const
{
    return new SyntaxForLoop( *this );
}


/** Get semantic value (not applicable so return NULL) */
RevPtr<RevVariable> SyntaxForLoop::evaluateContent( const std::shared_ptr<Environment>& env, bool dynamic )
{
    return NULL;
}


/** Finalize loop. */
void SyntaxForLoop::finalizeLoop( void )
{
    nextIndex = 0;
}


/** Get the name of the index variable */
const std::string& SyntaxForLoop::getIndexVarName( void ) const
{
    return varName;
}


/**
 * Get the next loop state. We simply assign a new value to
 * the loop variable using a control variable assignment.
 * We rely on the caller to make sure that the element index
 * is not out of range (or the container will throw an error).
 */
void SyntaxForLoop::getNextLoopState( void )
{
    // Get the next value from the container
    RevObject* elm = stateSpace->getElement( nextIndex-1 )->clone();

    // Set the loop variable to the next value using a control variable assignment
    loopVariable->replaceRevObject( elm );
    loopVariable->setWorkspaceVariableState( true );

    // Increment the element index
    nextIndex++;
}


/**
 * Initialize loop state. Here we evaluate the in-expression and
 * make sure it is a container with at least one dimension. We will
 * use the first dimension of the container to get the values for
 * the loop variable. We also add the loop variable to the environment.
 */
void SyntaxForLoop::initializeLoop( const std::shared_ptr<Environment>& env )
{
    
    assert ( nextIndex == 0 );  // Check that we are not running already

    // Evaluate expression and check that we get a vector
    const RevPtr<RevVariable>&      theVar   = inExpression->evaluateContent(env);
    const RevObject&             theValue    = theVar->getRevObject();

    // Check that it is a container (the first dimension of which we will use)
    stateSpace = dynamic_cast<Container*>( theValue.clone() );
    if ( stateSpace == NULL )
    {
       throw RbException( "The 'in' expression does not evaluate to a container" );
    }
    
    // Add the loop variable to the environment, if it is not already there
    if ( env->existsVariable( varName ) == false )
    {
        env->addVariable( varName, new RevVariable( NULL) );
    }
    
    // Set the local smart pointer to the loop variable
    loopVariable = env->getVariable( varName );
    
    // Initialize nextValue
    nextIndex = 1;
}


/**
 * Have we already reached the end of the loop? We return true
 * if we have reached the end.
 */
bool SyntaxForLoop::isFinished( void ) const
{
    return nextIndex > stateSpace->size();
}


/**
 * Is the syntax element safe for use in a function ( as
 * opposed to a procedure)? The variable element is safe
 * if the index variable is not an external variable, and
 * if the in-expression is function-safe.
 */
bool SyntaxForLoop::isFunctionSafe( const Environment& env, std::set<std::string>& localVars ) const
{
    if ( env.existsVariable( varName ) )
    {
        return false;
    }
    
    return inExpression->isFunctionSafe( env, localVars );
}

