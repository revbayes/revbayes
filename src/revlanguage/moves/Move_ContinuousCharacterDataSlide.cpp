#include <cstddef>
#include <ostream>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MetropolisHastingsMove.h"
#include "Move_ContinuousCharacterDataSlide.h"
#include "Natural.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlBoolean.h"
#include "RlContinuousCharacterData.h"
#include "TypeSpec.h"
#include "ContinuousCharacterDataSlideProposal.h"
#include "Move.h"
#include "RbBoolean.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlMove.h"
#include "StochasticNode.h"

namespace RevBayesCore { class ContinuousCharacterData; }
namespace RevBayesCore { class Proposal; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }


using namespace RevLanguage;

Move_ContinuousCharacterDataSlide::Move_ContinuousCharacterDataSlide() : Move()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_ContinuousCharacterDataSlide* Move_ContinuousCharacterDataSlide::clone(void) const
{
    
    return new Move_ContinuousCharacterDataSlide(*this);
}


void Move_ContinuousCharacterDataSlide::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new vector-slide move
    double l = static_cast<const RealPos &>( delta->getRevObject() ).getValue();
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    size_t del = static_cast<const Natural &>( delay->getRevObject() ).getValue();
    RevBayesCore::TypedDagNode<RevBayesCore::ContinuousCharacterData>* tmp = static_cast<const ContinuousCharacterData &>( x->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::ContinuousCharacterData> *n = static_cast<RevBayesCore::StochasticNode<RevBayesCore::ContinuousCharacterData> *>( tmp );
    bool t = static_cast<const RlBoolean &>( tune->getRevObject() ).getValue();
    
    RevBayesCore::Proposal *p = new RevBayesCore::ContinuousCharacterDataSlideProposal(n, l);
    value = new RevBayesCore::MetropolisHastingsMove(p, w, del, t);
    
}


/** Get Rev type of object */
const std::string& Move_ContinuousCharacterDataSlide::getClassType(void)
{
    
    static std::string rev_type = "Move_ContinuousCharacterDataSlide";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& Move_ContinuousCharacterDataSlide::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_ContinuousCharacterDataSlide::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "ContinuousCharacterDataSlide";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Move_ContinuousCharacterDataSlide::getParameterRules(void) const
{
    
    static MemberRules move_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        move_member_rules.push_back( new ArgumentRule( "x"     , ContinuousCharacterData::getClassTypeSpec(), "The variable this move operators on.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        move_member_rules.push_back( new ArgumentRule( "delta", RealPos::getClassTypeSpec()             , "The sliding window width parameter.", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new RealPos(1.0) ) );
        move_member_rules.push_back( new ArgumentRule( "tune"  , RlBoolean::getClassTypeSpec()           , "Should we tune this move during burnin?", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new RlBoolean( true ) ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return move_member_rules;
}

/** Get type spec */
const TypeSpec& Move_ContinuousCharacterDataSlide::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void Move_ContinuousCharacterDataSlide::printValue(std::ostream &o) const
{
    
    o << "Move_ContinuousCharacterDataSlide(";
    if (x != NULL)
    {
        o << x->getName();
    }
    else
    {
        o << "?";
    }
    o << ")";
}


/** Set a member variable */
void Move_ContinuousCharacterDataSlide::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "x" )
    {
        x = var;
    }
    else if ( name == "delta" )
    {
        delta = var;
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
