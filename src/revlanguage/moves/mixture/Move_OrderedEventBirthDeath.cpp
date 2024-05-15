#include <cstddef>
#include <ostream>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MetropolisHastingsMove.h"
#include "Move_OrderedEventBirthDeath.h"
#include "Natural.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlBoolean.h"
#include "RlString.h"
#include "TypeSpec.h"
#include "Move.h"
#include "OrderedEventTimes.h"
#include "OrderedEventBirthDeathProposal.h"
#include "RbBoolean.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlMove.h"
#include "RlOrderedEventTimes.h"
#include "StochasticNode.h"

namespace RevBayesCore { class OrderedEventTimes; }
namespace RevBayesCore { class Proposal; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }


using namespace RevLanguage;

/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Move_OrderedEventBirthDeath::Move_OrderedEventBirthDeath() : Move()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
Move_OrderedEventBirthDeath* Move_OrderedEventBirthDeath::clone(void) const
{
    
    return new Move_OrderedEventBirthDeath(*this);
}


/**
 * Create a new internal move object.
 *
 * This function simply dynamically allocates a new internal move object that is
 * associated with the variable (DAG-node). The internal move object is created by calling its
 * constructor and passing the move-parameters (the variable and other parameters) as arguments of the
 * constructor. The move constructor takes care of the proper hook-ups.
 *
 */
void Move_OrderedEventBirthDeath::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new random-geometric-walk move
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    size_t del = static_cast<const Natural &>( delay->getRevObject() ).getValue();
    RevBayesCore::TypedDagNode<RevBayesCore::OrderedEventTimes>* tmp = static_cast<const RlOrderedEventTimes &>( x->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::OrderedEventTimes> *n = static_cast<RevBayesCore::StochasticNode<RevBayesCore::OrderedEventTimes> *>( tmp );
    
    // finally create the internal move object
    RevBayesCore::Proposal *prop = new RevBayesCore::OrderedEventBirthDeathProposal(n);
    value = new RevBayesCore::MetropolisHastingsMove(prop, w, del, false);
}


/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Move_OrderedEventBirthDeath::getClassType(void)
{
    
    static std::string rev_type = "Move_OrderedEventBirthDeath";
    
    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Move_OrderedEventBirthDeath::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_OrderedEventBirthDeath::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "OrderedEventBirthDeath";
    
    return c_name;
}


/**
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the Slide move are:
 * (1) the variable which must be an integer.
 *
 * \return The member rules.
 */
const MemberRules& Move_OrderedEventBirthDeath::getParameterRules(void) const
{
    
    static MemberRules move_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        move_member_rules.push_back( new ArgumentRule( "x"     , RlOrderedEventTimes::getClassTypeSpec(),  "The variable on which this move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return move_member_rules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Move_OrderedEventBirthDeath::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}



void Move_OrderedEventBirthDeath::printValue(std::ostream &o) const
{
    
    o << "Move_OrderedEventBirthDeath(";
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


/**
 * Set a member variable.
 *
 * Sets a member variable with the given name and store the pointer to the variable.
 * The value of the variable might still change but this function needs to be called again if the pointer to
 * the variable changes. The current values will be used to create the distribution object.
 *
 * \param[in]    name     Name of the member variable.
 * \param[in]    var      Pointer to the variable.
 */
void Move_OrderedEventBirthDeath::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
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
