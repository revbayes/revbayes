#include <stddef.h>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MetropolisHastingsMove.h"
#include "ModelVector.h"
#include "Move_RandomCategoryWalk.h"
#include "Natural.h"
#include "RbException.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlBoolean.h"
#include "Probability.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "RandomCategoryWalkProposal.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "Move.h"
#include "RbBoolean.h"
#include "RbVector.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "RlMove.h"
#include "StochasticNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

namespace RevBayesCore { class Proposal; }


using namespace RevLanguage;

Move_RandomCategoryWalk::Move_RandomCategoryWalk() : Move()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_RandomCategoryWalk* Move_RandomCategoryWalk::clone(void) const
{
    
    return new Move_RandomCategoryWalk(*this);
}


void Move_RandomCategoryWalk::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // get the weight for this move
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    
    // get the stochastic variable
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<long> >* tmp = static_cast<const ModelVector<Natural> &>( x->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::RbVector<long> > *n = static_cast<RevBayesCore::StochasticNode<RevBayesCore::RbVector<long> > *>( tmp );
    
    // create the internal proposal object
    RevBayesCore::Proposal *prop = new RevBayesCore::RandomCategoryWalkProposal(n);
    
    // create the interval MH-move
    value = new RevBayesCore::MetropolisHastingsMove(prop,w,false);
}


/** Get Rev type of object */
const std::string& Move_RandomCategoryWalk::getClassType(void)
{
    
    static std::string rev_type = "Move_RandomCategoryWalk";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& Move_RandomCategoryWalk::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_RandomCategoryWalk::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "RandomCategoryWalk";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Move_RandomCategoryWalk::getParameterRules(void) const
{
    
    static MemberRules move_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        move_member_rules.push_back( new ArgumentRule( "x", ModelVector<Natural>::getClassTypeSpec()   , "The variable on which this move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return move_member_rules;
}

/** Get type spec */
const TypeSpec& Move_RandomCategoryWalk::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void Move_RandomCategoryWalk::printValue(std::ostream &o) const
{
    
    o << "Move_RandomCategoryWalk(";
    if (x != NULL)
    {
        o << x->getName();
    }
    else
    {
        o << "?";
    }
    o << ")";
}


/** Set a member variable */
void Move_RandomCategoryWalk::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
    if ( name == "x" )
    {
        x = var;
    }
    else
    {
        Move::setConstParameter(name, var);
    }
}
