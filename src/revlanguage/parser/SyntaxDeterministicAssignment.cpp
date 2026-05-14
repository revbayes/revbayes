#include <iostream>
#include <string>

#include "SyntaxDeterministicAssignment.h"
#include "DagNode.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "SyntaxAssignment.h"

namespace RevLanguage { class SyntaxElement; }

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
        // Name the rhs's DAG node (already in the graph) so cycle errors can show the variable.
        // The name is copied when the model clones the DAG.
        if ( !lhs->getName().empty() )
        {
            RevBayesCore::DagNode* rhs_node = rhs->getRevObject().getDagNode();
            if ( rhs_node != NULL )
            {
                rhs_node->setName( lhs->getName() );
            }
        }

        lhs->replaceRevObject( rhs->getRevObject().clone() );

        // Ensure the DAG node is named after the variable (e.g. "a" in "a := exp(b)") so cycle errors can identify it.
        RevBayesCore::DagNode* the_node = lhs->getRevObject().getDagNode();
        if ( the_node != NULL )
        {
            if ( !lhs->getName().empty() )
            {
                the_node->setName( lhs->getName() );
            }
            
            the_node->setParentNamePrefix( the_node->getName() );
        }
    }

}


/** Should we execute the rhs dynamically? Yes, because this is a deterministic assingment. */
bool SyntaxDeterministicAssignment::isDynamic( void )
{   
    return true;
}

