#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MetropolisHastingsMove.h"
#include "Move_GibbsBuddingTopology.h"
#include "GibbsBuddingTopologyProposal.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlTimeTree.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "Move.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlMove.h"
#include "StochasticNode.h"
#include "Tree.h"


using namespace RevLanguage;

Move_GibbsBuddingTopology::Move_GibbsBuddingTopology() : Move()
{

}


Move_GibbsBuddingTopology* Move_GibbsBuddingTopology::clone(void) const
{

    return new Move_GibbsBuddingTopology(*this);
}


void Move_GibbsBuddingTopology::constructInternalObject( void )
{
    // we free the memory first
    delete value;

    // now allocate a new move
    RevBayesCore::TypedDagNode<RevBayesCore::Tree> *tmp = static_cast<const TimeTree &>( tree->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::Tree> *t = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree> *>( tmp );

    double we = static_cast<const RealPos &>( weight->getRevObject() ).getValue();

    RevBayesCore::GibbsBuddingTopologyProposal *p = new RevBayesCore::GibbsBuddingTopologyProposal( t );

    value = new RevBayesCore::MetropolisHastingsMove(p, we, false);
}


/** Get Rev type of object */
const std::string& Move_GibbsBuddingTopology::getClassType(void)
{

    static std::string rev_type = "Move_GibbsBuddingTopology";

    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Move_GibbsBuddingTopology::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 */
std::string Move_GibbsBuddingTopology::getMoveName( void ) const
{
    std::string c_name = "GibbsBuddingTopology";

    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Move_GibbsBuddingTopology::getParameterRules(void) const
{

    static MemberRules move_member_rules;
    static bool rules_set = false;

    if ( rules_set == false )
    {
        move_member_rules.push_back( new ArgumentRule( "tree", TimeTree::getClassTypeSpec(), "The tree on which this move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );

        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );

        rules_set = true;
    }

    return move_member_rules;
}


/** Get type spec */
const TypeSpec& Move_GibbsBuddingTopology::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}


/** Get type spec */
void Move_GibbsBuddingTopology::printValue(std::ostream &o) const
{

    o << "Move_GibbsBuddingTopology(";
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


/** Set a member variable */
void Move_GibbsBuddingTopology::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
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
