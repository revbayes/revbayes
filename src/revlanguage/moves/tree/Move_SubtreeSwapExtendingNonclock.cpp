#include <stddef.h>
#include <ostream>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "SubtreeSwapExtending_nonClockProposal.h"
#include "RevObject.h"
#include "RealPos.h"
#include "MetropolisHastingsMove.h"
#include "Move_SubtreeSwapExtendingNonclock.h"
#include "RlBranchLengthTree.h"
#include "TypeSpec.h"
#include "Move.h"
#include "Probability.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlMove.h"
#include "StochasticNode.h"

namespace RevBayesCore { class Proposal; }
namespace RevBayesCore { class Tree; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }


using namespace RevLanguage;

Move_SubtreeSwapExtendingNonclock::Move_SubtreeSwapExtendingNonclock() : Move()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_SubtreeSwapExtendingNonclock* Move_SubtreeSwapExtendingNonclock::clone(void) const
{
    
	return new Move_SubtreeSwapExtendingNonclock(*this);
}


void Move_SubtreeSwapExtendingNonclock::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new sliding move
    double ep = static_cast<const RealPos &>( extension_prob->getRevObject() ).getValue();
    if ( ep >= 1.0 )
    {
        throw RbException( "The extensition probability should be strictly smaller than 1." );
    }
    
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    
    RevBayesCore::TypedDagNode<RevBayesCore::Tree> *tmp = static_cast<const BranchLengthTree &>( tree->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::Tree> *t = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree> *>( tmp );

    RevBayesCore::Proposal *p = new RevBayesCore::SubtreeSwapExtending_nonClockProposal(t, ep);
    value = new RevBayesCore::MetropolisHastingsMove(p, w);
    
}


/** Get class name of object */
const std::string& Move_SubtreeSwapExtendingNonclock::getClassName(void)
{
    
    static std::string rbClassName = "Move_SubtreeSwapExtending";
    
	return rbClassName;
}

/** Get class type spec describing type of object */
const TypeSpec& Move_SubtreeSwapExtendingNonclock::getClassTypeSpec(void)
{
    
    static TypeSpec rbClass = TypeSpec( getClassName(), new TypeSpec( Move::getClassTypeSpec() ) );
    
	return rbClass;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_SubtreeSwapExtendingNonclock::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "SubtreeSwapExtending";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Move_SubtreeSwapExtendingNonclock::getParameterRules(void) const
{
    
    static MemberRules SPRMemberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
    
        SPRMemberRules.push_back( new ArgumentRule( "tree", BranchLengthTree::getClassTypeSpec(), "The tree variable this move operates on.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        SPRMemberRules.push_back( new ArgumentRule( "extensionProbability", Probability::getClassTypeSpec(), "The extension probability (greater value, i.e., closer to 1, will result in bolder proposal).", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Probability(0.8) ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        SPRMemberRules.insert( SPRMemberRules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return SPRMemberRules;
}

/** Get type spec */
const TypeSpec& Move_SubtreeSwapExtendingNonclock::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}



/** Get type spec */
void Move_SubtreeSwapExtendingNonclock::printValue(std::ostream &o) const
{
    
    o << "SubtreeSwapExtending(";
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
void Move_SubtreeSwapExtendingNonclock::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "tree" )
    {
        tree = var;
    }
    else if ( name == "extensionProbability" )
    {
        extension_prob = var;
    }
    else
    {
        Move::setConstParameter(name, var);
    }
}
