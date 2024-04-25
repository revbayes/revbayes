#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "DirichletSimplexProposal.h"
#include "MetropolisHastingsMove.h"
#include "Move_DirichletSimplex.h"
#include "Natural.h"
#include "RevObject.h"
#include "Real.h"
#include "RealPos.h"
#include "RlBoolean.h"
#include "RlSimplex.h"
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

Move_DirichletSimplex::Move_DirichletSimplex() : Move()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_DirichletSimplex* Move_DirichletSimplex::clone(void) const
{
    
	return new Move_DirichletSimplex(*this);
}


void Move_DirichletSimplex::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new sliding move
    double a = static_cast<const RealPos &>( alpha->getRevObject() ).getValue();
    int nc   = static_cast<const Natural &>( numCats->getRevObject() ).getValue();
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    double o = static_cast<const RealPos &>( offset->getRevObject() ).getValue();
    double k = static_cast<const RealPos &>( kappa->getRevObject() ).getValue();
    double r = static_cast<const RealPos &>( tune_target->getRevObject() ).getValue();
    size_t del = static_cast<const Natural &>( delay->getRevObject() ).getValue();
    RevBayesCore::TypedDagNode< RevBayesCore::Simplex >* tmp = static_cast<const Simplex &>( x->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode< RevBayesCore::Simplex > *n = static_cast<RevBayesCore::StochasticNode< RevBayesCore::Simplex > *>( tmp );
    bool t = static_cast<const RlBoolean &>( tune->getRevObject() ).getValue();

    RevBayesCore::Proposal *prop = new RevBayesCore::DirichletSimplexProposal(n,a,nc,o,k,r);
    value = new RevBayesCore::MetropolisHastingsMove(prop,w,del,t);

}


/** Get Rev type of object */
const std::string& Move_DirichletSimplex::getClassType(void)
{
    
    static std::string rev_type = "Move_DirichletSimplex";
    
	return rev_type; 
}

/** Get class type spec describing type of object */
const TypeSpec& Move_DirichletSimplex::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Move_DirichletSimplex::getMoveAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "Simplex" );
    
    return a_names;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_DirichletSimplex::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "DirichletSimplex";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Move_DirichletSimplex::getParameterRules(void) const
{
    
    static MemberRules move_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        move_member_rules.push_back( new ArgumentRule( "x"      , Simplex::getClassTypeSpec()  , "The simplex on which this move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        move_member_rules.push_back( new ArgumentRule( "alpha"  , RealPos::getClassTypeSpec()  , "The concentration parameter on the previous value.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Real(1.0) ) );
        move_member_rules.push_back( new ArgumentRule( "numCats", Natural::getClassTypeSpec()  , "The number of categories changed per move.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(1) ) );
        move_member_rules.push_back( new ArgumentRule( "offset" , RealPos::getClassTypeSpec()  , "The offset of the current value to center new proposals (x+offset).", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RealPos(0.0) ) );
        move_member_rules.push_back( new ArgumentRule( "kappa"  , RealPos::getClassTypeSpec()  , "The offset of tuning parameter (x * alpha + kappa).", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RealPos(0.0) ) );
        move_member_rules.push_back( new ArgumentRule( "tune"   , RlBoolean::getClassTypeSpec(), "Should we tune this move during burnin?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( true ) ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return move_member_rules;
}

/** Get type spec */
const TypeSpec& Move_DirichletSimplex::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void Move_DirichletSimplex::printValue(std::ostream &o) const
{
    
    o << "Move_DirichletSimplex(";
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
void Move_DirichletSimplex::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
    if ( name == "x" )
    {
        x = var;
    }
    else if ( name == "alpha" )
    {
        alpha = var;
    }
    else if ( name == "numCats" )
    {
        numCats = var;
    }
    else if ( name == "offset" )
    {
        offset = var;
    }
    else if ( name == "kappa" )
    {
        kappa = var;
    }
    else if ( name == "tune" )
    {
        tune = var;
    }
    else {
        Move::setConstParameter(name, var);
    }
    
}
