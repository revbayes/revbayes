#include <cstddef>
#include <ostream>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RlBoolean.h"
#include "MetropolisHastingsMove.h"
#include "Move_ElementSwapSimplex.h"
#include "Natural.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlSimplex.h"
#include "ElementSwapSimplexProposal.h"
#include "TypeSpec.h"
#include "Move.h"
#include "RbBoolean.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlMove.h"
#include "StochasticNode.h"

namespace RevBayesCore { class Proposal; }
namespace RevBayesCore { class Simplex; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }


using namespace RevLanguage;

Move_ElementSwapSimplex::Move_ElementSwapSimplex() : Move()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_ElementSwapSimplex* Move_ElementSwapSimplex::clone(void) const
{
    
    return new Move_ElementSwapSimplex(*this);
}


void Move_ElementSwapSimplex::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new move
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    size_t del = static_cast<const Natural &>( delay->getRevObject() ).getValue();
    RevBayesCore::TypedDagNode< RevBayesCore::Simplex >* tmp = static_cast<const Simplex &>( x->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode< RevBayesCore::Simplex > *n = static_cast<RevBayesCore::StochasticNode< RevBayesCore::Simplex > *>( tmp );
    bool t = static_cast<const RlBoolean &>( tune->getRevObject() ).getValue();
    
    RevBayesCore::Proposal *prop = new RevBayesCore::ElementSwapSimplexProposal(n);
    value = new RevBayesCore::MetropolisHastingsMove(prop,w,del,t);
    
}


/** Get Rev type of object */
const std::string& Move_ElementSwapSimplex::getClassType(void)
{
    
    static std::string rev_type = "Move_ElementSwapSimplex";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& Move_ElementSwapSimplex::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_ElementSwapSimplex::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "ElementSwapSimplex";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Move_ElementSwapSimplex::getParameterRules(void) const
{
    
    static MemberRules move_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        move_member_rules.push_back( new ArgumentRule( "x"    , Simplex::getClassTypeSpec()  , "The variable this move operates on.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        move_member_rules.push_back( new ArgumentRule( "tune" , RlBoolean::getClassTypeSpec(), "Should we tune the concentration parameter during burnin?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( true ) ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return move_member_rules;
}

/** Get type spec */
const TypeSpec& Move_ElementSwapSimplex::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void Move_ElementSwapSimplex::printValue(std::ostream &o) const
{
    
    o << "Move_ElementSwapSimplex(";
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
void Move_ElementSwapSimplex::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "x" )
    {
        x = var;
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
