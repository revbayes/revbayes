#include <stddef.h>
#include <ostream>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "GibbsSubtreePruneAndRegraftProposal.h"
#include "RevObject.h"
#include "RealPos.h"
#include "MetropolisHastingsMove.h"
#include "Move_GibbsSubtreePruneAndRegraft.h"
#include "RlBranchLengthTree.h"
#include "TypeSpec.h"
#include "Move.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlMove.h"
#include "StochasticNode.h"

namespace RevBayesCore { class Proposal; }
namespace RevBayesCore { class Tree; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }


using namespace RevLanguage;

Move_GibbsSubtreePruneAndRegraft::Move_GibbsSubtreePruneAndRegraft() : Move()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_GibbsSubtreePruneAndRegraft* Move_GibbsSubtreePruneAndRegraft::clone(void) const
{
    
    return new Move_GibbsSubtreePruneAndRegraft(*this);
}


void Move_GibbsSubtreePruneAndRegraft::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new sliding move
    RevBayesCore::TypedDagNode<RevBayesCore::Tree> *tmp = static_cast<const BranchLengthTree &>( tree->getRevObject() ).getDagNode();
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    RevBayesCore::StochasticNode<RevBayesCore::Tree> *t = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree> *>( tmp );

    RevBayesCore::Proposal *p = new RevBayesCore::GibbsSubtreePruneAndRegraftProposal(t);
    value = new RevBayesCore::MetropolisHastingsMove(p,w);
    
}


/** Get class name of object */
const std::string& Move_GibbsSubtreePruneAndRegraft::getClassName(void)
{
    
    static std::string rbClassName = "Move_GibbsSubtreePruneAndRegraft";
    
    return rbClassName;
}

/** Get class type spec describing type of object */
const TypeSpec& Move_GibbsSubtreePruneAndRegraft::getClassTypeSpec(void)
{
    
    static TypeSpec rbClass = TypeSpec( getClassName(), new TypeSpec( Move::getClassTypeSpec() ) );
    
    return rbClass;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_GibbsSubtreePruneAndRegraft::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "GibbsSubtreePruneAndRegraft";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Move_GibbsSubtreePruneAndRegraft::getParameterRules(void) const
{
    
    static MemberRules move_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
    
        move_member_rules.push_back( new ArgumentRule( "tree", BranchLengthTree::getClassTypeSpec(), "The tree variable this move operates on.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return move_member_rules;
}

/** Get type spec */
const TypeSpec& Move_GibbsSubtreePruneAndRegraft::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}



/** Get type spec */
void Move_GibbsSubtreePruneAndRegraft::printValue(std::ostream &o) const
{
    
    o << "GibbsSubtreePruneAndRegraft(";
    if (tree != NULL)
    {
        o << tree->getName();
    }
    else
    {
        o << "?";
    }
    o << ")";
}


/** Set a NearestNeighborInterchange variable */
void Move_GibbsSubtreePruneAndRegraft::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "tree" )
    {
        tree = var;
    }
    else
    {
        Move::setConstParameter(name, var);
    }
}
