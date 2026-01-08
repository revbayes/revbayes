#include <cstddef>
#include <ostream>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MetropolisHastingsMove.h"
#include "Move_RootTimeScaleBactrian.h"
#include "Natural.h"
#include "RootTimeScaleBactrianProposal.h"
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

Move_RootTimeScaleBactrian::Move_RootTimeScaleBactrian() : Move()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_RootTimeScaleBactrian* Move_RootTimeScaleBactrian::clone(void) const
{
    
    return new Move_RootTimeScaleBactrian(*this);
}


void Move_RootTimeScaleBactrian::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new move
    RevBayesCore::TypedDagNode<RevBayesCore::Tree> *tmp = static_cast<const TimeTree &>( tree->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::Tree> *t = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree> *>( tmp );
    
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    size_t del = static_cast<const Natural &>( delay->getRevObject() ).getValue();
    
    RevBayesCore::Proposal *p = new RevBayesCore::RootTimeScaleBactrianProposal( t );
    value = new RevBayesCore::MetropolisHastingsMove(p,w,del,false);
}


/** Get Rev type of object */
const std::string& Move_RootTimeScaleBactrian::getClassType(void)
{
    
    static std::string rev_type = "Move_RootTimeScaleBactrian";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& Move_RootTimeScaleBactrian::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_RootTimeScaleBactrian::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "RootTimeScaleBactrian";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Move_RootTimeScaleBactrian::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        memberRules.push_back( new ArgumentRule( "tree", TimeTree::getClassTypeSpec(), "The tree on which this move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        memberRules.insert( memberRules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return memberRules;
}

/** Get type spec */
const TypeSpec& Move_RootTimeScaleBactrian::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}



/** Get type spec */
void Move_RootTimeScaleBactrian::printValue(std::ostream &o) const
{
    
    o << "Move_RootTimeScaleBactrian(";
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
void Move_RootTimeScaleBactrian::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
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
