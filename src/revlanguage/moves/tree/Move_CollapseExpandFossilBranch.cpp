#include <cstddef>
#include <ostream>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MetropolisHastingsMove.h"
#include "Move_CollapseExpandFossilBranch.h"
#include "CollapseExpandFossilBranchProposal.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlTimeTree.h"
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

Move_CollapseExpandFossilBranch::Move_CollapseExpandFossilBranch() : Move()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_CollapseExpandFossilBranch* Move_CollapseExpandFossilBranch::clone(void) const
{
    
    return new Move_CollapseExpandFossilBranch(*this);
}


void Move_CollapseExpandFossilBranch::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new move
    RevBayesCore::TypedDagNode<RevBayesCore::Tree> *tmp = static_cast<const TimeTree &>( tree->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::Tree> *t = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree> *>( tmp );
    
    RevBayesCore::TypedDagNode<double> *o = static_cast<const RealPos &>( origin->getRevObject() ).getDagNode();
    
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    
    RevBayesCore::Proposal *p = new RevBayesCore::CollapseExpandFossilBranchProposal( t, o );
    value = new RevBayesCore::MetropolisHastingsMove(p,w,false);
}


/** Get Rev type of object */
const std::string& Move_CollapseExpandFossilBranch::getClassType(void)
{
    
    static std::string rev_type = "Move_CollapseExpandFossilBranch";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& Move_CollapseExpandFossilBranch::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_CollapseExpandFossilBranch::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "CollapseExpandFossilBranch";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Move_CollapseExpandFossilBranch::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        memberRules.push_back( new ArgumentRule( "tree"  , TimeTree::getClassTypeSpec(), "The tree on which this moves operates on. It should be a fossil tree!", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        memberRules.push_back( new ArgumentRule( "origin", RealPos::getClassTypeSpec() , "The variable for the origin of the process giving a maximum age for the new fossil attachement time.", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        memberRules.insert( memberRules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return memberRules;
}

/** Get type spec */
const TypeSpec& Move_CollapseExpandFossilBranch::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}



/** Get type spec */
void Move_CollapseExpandFossilBranch::printValue(std::ostream &o) const
{
    
    o << "Move_CollapseExpandFossilBranch(";
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
void Move_CollapseExpandFossilBranch::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "tree" )
    {
        tree = var;
    }
    else if ( name == "origin" )
    {
        origin = var;
    }
    else
    {
        Move::setConstParameter(name, var);
    }
}
