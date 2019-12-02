#include "RbException.h"
#include "RbOptions.h"
#include "SyntaxDeterministicAssignment.h"
#include "Workspace.h"

#include <iostream>
#include <list>
#include <sstream>

using namespace RevLanguage;

/** Construct from operator type, variable and expression */
SyntaxDeterministicAssignment::SyntaxDeterministicAssignment( SyntaxElement* lhsExpr, SyntaxElement* rhsExpr ) : SyntaxAssignment(lhsExpr,rhsExpr)
{

}


/** Destructor deletes operands */
SyntaxDeterministicAssignment::~SyntaxDeterministicAssignment()
{

}


/** Type-safe clone of syntax element */
SyntaxDeterministicAssignment* SyntaxDeterministicAssignment::clone () const
{
    return new SyntaxDeterministicAssignment( *this );
}


/** Get semantic value: insert symbol and return the rhs value of the assignment */
void SyntaxDeterministicAssignment::assign(RevPtr<RevVariable> &lhs, RevPtr<RevVariable> &rhs)
{

    // Check if the variable returned from the rhs expression is a named
    // variable in the environment. If so, we want to create an indirect
    // reference to it; otherwise, we want to fill the slot with a clone
    // of the variable returned by the rhs expression.
    if ( rhs->getName() != "" )
    {
        lhs->replaceRevObject( rhs->getRevObject().makeIndirectReference() );
    }
    else
    {    
        lhs->replaceRevObject( rhs->getRevObject().clone() );
        
        // make sure all the implicitly created variables got a correct name
        RevBayesCore::DagNode* the_node = lhs->getRevObject().getDagNode();
        the_node->setParentNamePrefix( the_node->getName() );
    }

}


/** Should we execute the rhs dynamically? Yes, because this is a deterministic assingment. */
bool SyntaxDeterministicAssignment::isDynamic( void )
{   
    return true;
}

