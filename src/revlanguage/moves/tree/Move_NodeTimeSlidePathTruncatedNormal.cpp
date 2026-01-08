#include <cstddef>
#include <ostream>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MetropolisHastingsMove.h"
#include "Move_NodeTimeSlidePathTruncatedNormal.h"
#include "Natural.h"
#include "NodeTimeSlidePathTruncatedNormalProposal.h"
#include "RlBoolean.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlTimeTree.h"
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

Move_NodeTimeSlidePathTruncatedNormal::Move_NodeTimeSlidePathTruncatedNormal() : Move()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_NodeTimeSlidePathTruncatedNormal* Move_NodeTimeSlidePathTruncatedNormal::clone(void) const
{
    
	return new Move_NodeTimeSlidePathTruncatedNormal(*this);
}


void Move_NodeTimeSlidePathTruncatedNormal::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new move
    RevBayesCore::TypedDagNode<RevBayesCore::Tree> *tmp = static_cast<const TimeTree &>( tree->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::Tree> *tau = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree> *>( tmp );

    double s = static_cast<const RealPos &>( sigma->getRevObject() ).getValue();
    
    bool a = static_cast<const RlBoolean &>( scaleByAge->getRevObject() ).getValue();

    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    
    size_t del = static_cast<const Natural &>( delay->getRevObject() ).getValue();
    
    bool t = static_cast<const RlBoolean &>( tune->getRevObject() ).getValue();

    RevBayesCore::Proposal *p = new RevBayesCore::NodeTimeSlidePathTruncatedNormalProposal( tau, s, a );
    value = new RevBayesCore::MetropolisHastingsMove(p,w,del,t);
}


/** Get Rev type of object */
const std::string& Move_NodeTimeSlidePathTruncatedNormal::getClassType(void)
{
    
    static std::string rev_type = "Move_NodeTimeSlidePathTruncatedNormal";
    
	return rev_type; 
}

/** Get class type spec describing type of object */
const TypeSpec& Move_NodeTimeSlidePathTruncatedNormal::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_NodeTimeSlidePathTruncatedNormal::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "NodeTimeSlidePathTruncatedNormal";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Move_NodeTimeSlidePathTruncatedNormal::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        memberRules.push_back( new ArgumentRule( "tree", TimeTree::getClassTypeSpec(), "The tree on which this move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        memberRules.push_back( new ArgumentRule( "sigma", TimeTree::getClassTypeSpec(), "The sd of the proposal.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC, new RealPos( 1.0 )) );
        memberRules.push_back( new ArgumentRule( "scaleByAge"  , RlBoolean::getClassTypeSpec(), "Should we additionally use age of node to adjust sd?", ArgumentRule::BY_VALUE    , ArgumentRule::ANY       , new RlBoolean( true ) ) );
        memberRules.push_back( new ArgumentRule( "tune"  , RlBoolean::getClassTypeSpec(), "Should we tune sigma during burnin?", ArgumentRule::BY_VALUE    , ArgumentRule::ANY       , new RlBoolean( true ) ) );

        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        memberRules.insert( memberRules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return memberRules;
}

/** Get type spec */
const TypeSpec& Move_NodeTimeSlidePathTruncatedNormal::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}



/** Get type spec */
void Move_NodeTimeSlidePathTruncatedNormal::printValue(std::ostream &o) const
{
    
    o << "Move_NodeTimeSlidePathTruncatedNormal(";
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
void Move_NodeTimeSlidePathTruncatedNormal::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{

    if ( name == "tree" )
    {
        tree = var;
    }
    else if ( name == "sigma" )
    {
        sigma = var;
    }
    else if ( name == "scaleByAge" )
    {
        scaleByAge = var;
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
