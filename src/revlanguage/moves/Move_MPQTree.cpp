#include "Move_MPQTree.h"

#include <cstddef>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MPQRateMatrixProposal.h"
#include "RlBoolean.h"
#include "MetropolisHastingsMove.h"
#include "ModelVector.h"
#include "Real.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlRateGenerator.h"
#include "RlTree.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "ModelObject.h"
#include "Move.h"
#include "RbBoolean.h"
#include "StochasticNode.h"
#include "StringUtilities.h"

namespace RevBayesCore { class MatrixReal; }
namespace RevBayesCore { class Proposal; }
namespace RevBayesCore { template <class valueType> class RbVector; }


using namespace RevLanguage;

Move_MPQTree::Move_MPQTree() : Move() {
    
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_MPQTree* Move_MPQTree::clone(void) const {
    
    return new Move_MPQTree(*this);
}


void Move_MPQTree::constructInternalObject( void ) {

    // we free the memory first
    delete value;
    
    // now allocate a new sliding move
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    RevBayesCore::TypedDagNode<RevBayesCore::RateGenerator >* tmp = static_cast<const RateGenerator &>( Q->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::RateGenerator > *n = static_cast<RevBayesCore::StochasticNode<RevBayesCore::RateGenerator> *>( tmp );

    RevBayesCore::TypedDagNode<RevBayesCore::Tree >* tmp_tree = static_cast<const Tree &>( tree->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::Tree > *t = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree> *>( tmp_tree );

    bool tu = static_cast<const RlBoolean &>( tune->getRevObject() ).getValue();
    
    RevBayesCore::Proposal *p = new RevBayesCore::MPQRateMatrixProposal(n);
    value = new RevBayesCore::MetropolisHastingsMove(p,w,tu);

}


/** Get class name of object */
const std::string& Move_MPQTree::getClassType(void) {
    
    static std::string revClassType = "Move_MPQTree";
    
    return revClassType;
}

/** Get class type spec describing type of object */
const TypeSpec& Move_MPQTree::getClassTypeSpec(void) {
    
    static TypeSpec revClassTypeSpec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
    return revClassTypeSpec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_MPQTree::getMoveName( void ) const {

    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "MPQRateMatrix";
    
    return c_name;
}



/** Return member rules (no members) */
const MemberRules& Move_MPQTree::getParameterRules(void) const {
    
    static MemberRules move_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
        {
        move_member_rules.push_back( new ArgumentRule( "Q"     , RateGenerator::getClassTypeSpec(),        "The general nucleotide rate matrix on which this move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        move_member_rules.push_back( new ArgumentRule( "tree"  , Tree::getClassTypeSpec(),                 "The phylogeny.", ArgumentRule::BY_VALUE    , ArgumentRule::ANY ) );
        move_member_rules.push_back( new ArgumentRule( "tune"  , RlBoolean::getClassTypeSpec(),            "Should we tune the scaling factor during burnin?", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new RlBoolean( true ) ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
        }
    
    return move_member_rules;
}

/** Get type spec */
const TypeSpec& Move_MPQTree::getTypeSpec( void ) const {
    
    static TypeSpec type_spec = getClassTypeSpec();
    return type_spec;
}


/** Get type spec */
void Move_MPQTree::printValue(std::ostream &o) const {
    
    o << "Move_MPQTree(";
    if (Q != NULL)
        {
        o << Q->getName();
        }
    else
        {
        o << "?";
        }
    o << ")";
}


/** Set a member variable */
void Move_MPQTree::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
    if ( name == "Q" )
        {
        Q = var;
        }
    else if ( name == "tree" )
        {
        tree = var;
        }
    else if ( name == "weight" )
        {
        weight = var;
        }
    else if ( name == "tune" )
        {
        tune = var;
        }
    else
        {
        Move::setConstParameter(name, var);
        }
}

