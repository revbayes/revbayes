#include <stddef.h>
#include <ostream>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Natural.h"
#include "MetropolisHastingsMove.h"
#include "Move_WeightedBranchLengthScale.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlBoolean.h"
#include "RlBranchLengthTree.h"
#include "WeightedBranchLengthScaleProposal.h"
#include "TypeSpec.h"
#include "Move.h"
#include "RbBoolean.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlMove.h"
#include "StochasticNode.h"

namespace RevBayesCore { class Proposal; }
namespace RevBayesCore { class Tree; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }


using namespace RevLanguage;

Move_WeightedBranchLengthScale::Move_WeightedBranchLengthScale() : Move()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_WeightedBranchLengthScale* Move_WeightedBranchLengthScale::clone(void) const
{
    
    return new Move_WeightedBranchLengthScale(*this);
}


void Move_WeightedBranchLengthScale::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new tree scale move
    RevBayesCore::TypedDagNode<RevBayesCore::Tree> *tmp = static_cast<const BranchLengthTree &>( tree->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::Tree> *t = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree> *>( tmp );
    
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    long n = static_cast<const Natural &>( num_breaks->getRevObject() ).getValue();
    double a = static_cast<const RealPos &>( alpha->getRevObject() ).getValue();
//    bool tune = static_cast<const RlBoolean &>( tuning->getRevObject() ).getValue();
    bool tune = false;
    
    RevBayesCore::Proposal *p = new RevBayesCore::WeightedBranchLengthScaleProposal(t, n, a);
    value = new RevBayesCore::MetropolisHastingsMove(p, w, tune);
}


/** Get Rev type of object */
const std::string& Move_WeightedBranchLengthScale::getClassType(void)
{
    
    static std::string rev_type = "Move_WeightedBranchLengthScale";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& Move_WeightedBranchLengthScale::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_WeightedBranchLengthScale::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "WeightedBranchLengthScale";
    
    return c_name;
}


/** Return member rules */
const MemberRules& Move_WeightedBranchLengthScale::getParameterRules(void) const
{
    
    static MemberRules move_member_rules;
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        move_member_rules.push_back( new ArgumentRule( "tree"       , BranchLengthTree::getClassTypeSpec() , "The tree variable the move operates on.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        move_member_rules.push_back( new ArgumentRule( "alpha"      , RealPos::getClassTypeSpec()  , "The parameter of the beta distribution to compute the break points.", ArgumentRule::BY_VALUE    , ArgumentRule::ANY       , new RealPos( 0.5 ) ) );
        move_member_rules.push_back( new ArgumentRule( "numBreaks"  , Natural::getClassTypeSpec(), "The number of break points", ArgumentRule::BY_VALUE    , ArgumentRule::ANY       , new Natural( 6 ) ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return move_member_rules;
}

/** Get type spec */
const TypeSpec& Move_WeightedBranchLengthScale::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}



/** Get type spec */
void Move_WeightedBranchLengthScale::printValue(std::ostream &o) const
{
    
    o << "Move_WeightedBranchLengthScale(";
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
void Move_WeightedBranchLengthScale::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "tree" )
    {
        tree = var;
    }
    else if ( name == "alpha" )
    {
        alpha = var;
    }
    else if ( name == "numBreaks" )
    {
        num_breaks = var;
    }
    else
    {
        Move::setConstParameter(name, var);
    }
    
}

