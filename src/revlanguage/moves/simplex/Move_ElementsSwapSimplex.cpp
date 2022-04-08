#include <stddef.h>
#include <ostream>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RlBoolean.h"
#include "MetropolisHastingsMove.h"
#include "Move_ElementsSwapSimplex.h"
#include "Natural.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlSimplex.h"
#include "ElementsSwapSimplexProposal.h"
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

Move_ElementsSwapSimplex::Move_ElementsSwapSimplex() : Move()
{
    // first, the argument rules
    ArgumentRules* addIndexRules = new ArgumentRules();
    
    // next, set the specific arguments
    addIndexRules->push_back( new ArgumentRule( "elementIndex", Natural::getClassTypeSpec(), "The index of the element to add to the swapping set.", ArgumentRule::BY_VALUE , ArgumentRule::ANY, new Natural( 1 ) ) );
    
    // finally, create the methods
    methods.addFunction( new MemberProcedure( "addIndex", RlUtils::Void, addIndexRules) );
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_ElementsSwapSimplex* Move_ElementsSwapSimplex::clone(void) const
{
    
    return new Move_ElementsSwapSimplex(*this);
}


void Move_ElementsSwapSimplex::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new move
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    RevBayesCore::TypedDagNode< RevBayesCore::Simplex >* tmp = static_cast<const Simplex &>( x->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode< RevBayesCore::Simplex > *n = static_cast<RevBayesCore::StochasticNode< RevBayesCore::Simplex > *>( tmp );
    
    RevBayesCore::Proposal *prop = new RevBayesCore::ElementsSwapSimplexProposal(n);
    value = new RevBayesCore::MetropolisHastingsMove(prop, w);
    
}


RevPtr<RevVariable> Move_ElementsSwapSimplex::executeMethod(const std::string& name, const std::vector<Argument>& args, bool &found)
{
    
    if ( name == "addIndex" )
    {
        found = true;
        
        size_t idx = static_cast<const Natural &>( args[0].getVariable()->getRevObject() ).getValue();
        
        RevBayesCore::MetropolisHastingsMove *m = static_cast<RevBayesCore::MetropolisHastingsMove*>(this->value);
        RevBayesCore::ElementsSwapSimplexProposal &prop = static_cast<RevBayesCore::ElementsSwapSimplexProposal&>( m->getProposal() );
        
        RevBayesCore::TypedDagNode< RevBayesCore::Simplex>* tmp = static_cast<const Simplex &>( x->getRevObject() ).getDagNode();
        RevBayesCore::StochasticNode< RevBayesCore::Simplex> *n = static_cast<RevBayesCore::StochasticNode< RevBayesCore::Simplex > *>( tmp );
        
        if ( idx >= 1 && idx <= n->getValue().size())
        {
            prop.addIndex( idx - 1 );
        }
        else
        {
            std::stringstream ss_err;
            ss_err << "Index out of bounds: The vector of size " << n->getValue().size() << " does not have an element for index " << idx << ".";
            throw RbException(ss_err.str());
        }
        
        return NULL;
    }
    
    return Move::executeMethod( name, args, found );
}


/** Get Rev type of object */
const std::string& Move_ElementsSwapSimplex::getClassType(void)
{
    
    static std::string rev_type = "Move_ElementsSwapSimplex";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Move_ElementsSwapSimplex::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_ElementsSwapSimplex::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "ElementsSwapSimplex";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Move_ElementsSwapSimplex::getParameterRules(void) const
{
    
    static MemberRules move_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        move_member_rules.push_back( new ArgumentRule( "x"    , Simplex::getClassTypeSpec()  , "The variable this move operates on.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return move_member_rules;
}

/** Get type spec */
const TypeSpec& Move_ElementsSwapSimplex::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void Move_ElementsSwapSimplex::printValue(std::ostream &o) const
{
    
    o << "Move_ElementsSwapSimplex(";
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
void Move_ElementsSwapSimplex::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "x" )
    {
        x = var;
    }
    else
    {
        Move::setConstParameter(name, var);
    }
}
