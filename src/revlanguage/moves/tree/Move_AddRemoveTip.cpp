#include <cstddef>
#include <ostream>
#include <string>

#include "AddRemoveTipProposal.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MetropolisHastingsMove.h"
#include "Move_AddRemoveTip.h"
#include "RbException.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlBoolean.h"
#include "RlTimeTree.h"
#include "TypeSpec.h"
#include "Move.h"
#include "RbBoolean.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlMove.h"
#include "StochasticNode.h"
#include "StringUtilities.h"

namespace RevBayesCore { class Proposal; }
namespace RevBayesCore { class Tree; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }


using namespace RevLanguage;

Move_AddRemoveTip::Move_AddRemoveTip() : Move()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_AddRemoveTip* Move_AddRemoveTip::clone(void) const
{
    
    return new Move_AddRemoveTip(*this);
}


void Move_AddRemoveTip::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new move
    RevBayesCore::TypedDagNode<RevBayesCore::Tree> *tmp = static_cast<const TimeTree &>( tree->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::Tree> *t = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree> *>( tmp );
    
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    
    bool exa = static_cast<const RlBoolean &>( extant->getRevObject() ).getValue();

    bool exi = static_cast<const RlBoolean &>( extinct->getRevObject() ).getValue();

    bool sa = static_cast<const RlBoolean &>( sampled_ancestors->getRevObject() ).getValue();

    if (exa == false && exi == false)
    {
        throw(RbException("In mvAddRemoveTip, 'extant' and 'extinct' cannot both be false"));
    }

    RevBayesCore::Proposal *p = new RevBayesCore::AddRemoveTipProposal( t, exa, exi, sa );
    value = new RevBayesCore::MetropolisHastingsMove(p,w,false);
}


/** Get Rev type of object */
const std::string& Move_AddRemoveTip::getClassType(void)
{
    
    static std::string rev_type = "Move_AddRemoveTip";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& Move_AddRemoveTip::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_AddRemoveTip::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "AddRemoveTip";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Move_AddRemoveTip::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        memberRules.push_back( new ArgumentRule( "tree"  , TimeTree::getClassTypeSpec(), "The tree on which this moves operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        memberRules.push_back( new ArgumentRule( "extant", RlBoolean::getClassTypeSpec(), "Should we add/remove extant tips?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( true ) ) );
        memberRules.push_back( new ArgumentRule( "extinct", RlBoolean::getClassTypeSpec(), "Should we add/remove extinct tips?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( false ) ) );
        memberRules.push_back( new ArgumentRule( "sa", RlBoolean::getClassTypeSpec(), "Should we add/remove extinct tips as sampled ancestors?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( false ) ) );


        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        memberRules.insert( memberRules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return memberRules;
}

/** Get type spec */
const TypeSpec& Move_AddRemoveTip::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}



/** Get type spec */
void Move_AddRemoveTip::printValue(std::ostream &o) const
{
    
    o << "Move_AddRemoveTip(";
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
void Move_AddRemoveTip::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "tree" )
    {
        tree = var;
    }
    else if ( name == "extant" )
    {
        extant = var;
    }
    else if ( name == "extinct" )
    {
        extinct = var;
    }
    else if ( name == "sa" )
    {
        sampled_ancestors = var;
    }
    else
    {
        Move::setConstParameter(name, var);
    }
}
