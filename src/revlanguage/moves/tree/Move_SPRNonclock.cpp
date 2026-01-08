#include <cstddef>
#include <ostream>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "SubtreePruneRegraftProposal.h"
#include "RevObject.h"
#include "RealPos.h"
#include "MetropolisHastingsMove.h"
#include "Move_SPRNonclock.h"
#include "Natural.h"
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

Move_SPRNonclock::Move_SPRNonclock() : Move()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_SPRNonclock* Move_SPRNonclock::clone(void) const
{
    
	return new Move_SPRNonclock(*this);
}


void Move_SPRNonclock::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new sliding move
    RevBayesCore::TypedDagNode<RevBayesCore::Tree> *tmp = static_cast<const BranchLengthTree &>( tree->getRevObject() ).getDagNode();
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    size_t del = static_cast<const Natural &>( delay->getRevObject() ).getValue();
    RevBayesCore::StochasticNode<RevBayesCore::Tree> *t = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree> *>( tmp );

    RevBayesCore::Proposal *p = new RevBayesCore::SubtreePruneRegraftProposal(t);
    value = new RevBayesCore::MetropolisHastingsMove(p,w,del,false);
    
}


/** Get class name of object */
const std::string& Move_SPRNonclock::getClassName(void)
{
    
    static std::string rbClassName = "Move_SPR";
    
	return rbClassName;
}

/** Get class type spec describing type of object */
const TypeSpec& Move_SPRNonclock::getClassTypeSpec(void)
{
    
    static TypeSpec rbClass = TypeSpec( getClassName(), new TypeSpec( Move::getClassTypeSpec() ) );
    
	return rbClass;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_SPRNonclock::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "SPR";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Move_SPRNonclock::getParameterRules(void) const
{
    
    static MemberRules SPRMemberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
    
        SPRMemberRules.push_back( new ArgumentRule( "tree", BranchLengthTree::getClassTypeSpec(), "The tree variable this move operates on.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        SPRMemberRules.insert( SPRMemberRules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return SPRMemberRules;
}

/** Get type spec */
const TypeSpec& Move_SPRNonclock::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}



/** Get type spec */
void Move_SPRNonclock::printValue(std::ostream &o) const
{
    
    o << "SPR(";
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
void Move_SPRNonclock::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
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
