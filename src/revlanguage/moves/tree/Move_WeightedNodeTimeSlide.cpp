#include <cstddef>
#include <ostream>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MetropolisHastingsMove.h"
#include "Move_WeightedNodeTimeSlide.h"
#include "Natural.h"
#include "NodeTimeSlideWeightedProposal.h"
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

Move_WeightedNodeTimeSlide::Move_WeightedNodeTimeSlide() : Move()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_WeightedNodeTimeSlide* Move_WeightedNodeTimeSlide::clone(void) const
{
    
	return new Move_WeightedNodeTimeSlide(*this);
}


void Move_WeightedNodeTimeSlide::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new sliding move
    RevBayesCore::TypedDagNode<RevBayesCore::Tree> *tmp = static_cast<const TimeTree &>( tree->getRevObject() ).getDagNode();
    size_t b = static_cast<const Natural &>( blocks->getRevObject() ).getValue();
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    size_t del = static_cast<const Natural &>( delay->getRevObject() ).getValue();
    RevBayesCore::StochasticNode<RevBayesCore::Tree> *t = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree> *>( tmp );
    
    RevBayesCore::Proposal *p = new RevBayesCore::NodeTimeSlideWeightedProposal(t, b);
    value = new RevBayesCore::MetropolisHastingsMove(p, w, del, false);

}


/** Get Rev type of object */
const std::string& Move_WeightedNodeTimeSlide::getClassType(void) { 
    
    static std::string rev_type = "Move_WeightedNodeTimeSlide";
    
	return rev_type; 
}

/** Get class type spec describing type of object */
const TypeSpec& Move_WeightedNodeTimeSlide::getClassTypeSpec(void) { 
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_WeightedNodeTimeSlide::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "NodeTimeSlide";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Move_WeightedNodeTimeSlide::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        memberRules.push_back( new ArgumentRule( "tree"  , TimeTree::getClassTypeSpec(), "The tree on which this move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        memberRules.push_back( new ArgumentRule( "blocks", Natural::getClassTypeSpec() , "The number of bocks into which the branch will be broken.", ArgumentRule::BY_VALUE    , ArgumentRule::ANY       , new Natural(8) ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        memberRules.insert( memberRules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return memberRules;
}

/** Get type spec */
const TypeSpec& Move_WeightedNodeTimeSlide::getTypeSpec( void ) const {
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}



/** Get type spec */
void Move_WeightedNodeTimeSlide::printValue(std::ostream &o) const {
    
    o << "Move_WeightedNodeTimeSlide(";
    if (tree != NULL) {
        o << tree->getName();
    }
    else {
        o << "?";
    }
    o << ")";
}


/** Set a NearestNeighborInterchange variable */
void Move_WeightedNodeTimeSlide::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
    if ( name == "tree" ) {
        tree = var;
    }
    else if ( name == "blocks" ) {
        blocks = var;
    }
    else {
        Move::setConstParameter(name, var);
    }
}
