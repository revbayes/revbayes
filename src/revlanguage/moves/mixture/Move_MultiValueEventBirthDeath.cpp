#include <cstddef>
#include <ostream>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MetropolisHastingsMove.h"
#include "Move_MultiValueEventBirthDeath.h"
#include "MultiValueEventBirthDeathProposal.h"
#include "Natural.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlBoolean.h"
#include "TypeSpec.h"
#include "Move.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlMove.h"
#include "RlMultiValueEvent.h"
#include "StochasticNode.h"

namespace RevBayesCore { class MultiValueEvent; }
namespace RevBayesCore { class Proposal; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }


using namespace RevLanguage;

/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Move_MultiValueEventBirthDeath::Move_MultiValueEventBirthDeath() : Move()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
Move_MultiValueEventBirthDeath* Move_MultiValueEventBirthDeath::clone(void) const
{
    
    return new Move_MultiValueEventBirthDeath(*this);
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
void Move_MultiValueEventBirthDeath::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new random-geometric-walk move
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    size_t del = static_cast<const Natural &>( delay->getRevObject() ).getValue();
    RevBayesCore::TypedDagNode<RevBayesCore::MultiValueEvent>* tmp = static_cast<const MultiValueEvent &>( x->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::MultiValueEvent> *n = static_cast<RevBayesCore::StochasticNode<RevBayesCore::MultiValueEvent> *>( tmp );
    bool use_ac = static_cast<const RlBoolean &>( ac->getRevObject() ).getValue();
//    bool t = static_cast<const RlBoolean &>( tune->getRevObject() ).getValue();
    bool t = false;

    // finally create the internal move object
    RevBayesCore::Proposal *prop = new RevBayesCore::MultiValueEventBirthDeathProposal(n,use_ac);
    value = new RevBayesCore::MetropolisHastingsMove(prop,w,del,t);
    
}


/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Move_MultiValueEventBirthDeath::getClassType(void)
{
    
    static std::string rev_type = "Move_MultiValueEventBirthDeath";
    
    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Move_MultiValueEventBirthDeath::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_MultiValueEventBirthDeath::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "MultiValueEventBirthDeath";
    
    return c_name;
}


/**
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the scale move are:
 * (1) the variable which must be an integer.
 *
 * \return The member rules.
 */
const MemberRules& Move_MultiValueEventBirthDeath::getParameterRules(void) const
{
    
    static MemberRules move_member_rules;
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        
        move_member_rules.push_back( new ArgumentRule( "x"   , MultiValueEvent::getClassTypeSpec(),  "The variable on which this move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        move_member_rules.push_back( new ArgumentRule( "ac", RlBoolean::getClassTypeSpec(), "Should we use the autocorrelated proposal?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( false ) ) );
//        move_member_rules.push_back( new ArgumentRule( "tune", RlBoolean::getClassTypeSpec(), "Should we tune the scaling factor during burnin?", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new RlBoolean( true ) ) );

        /* Inherit weight from Move, put it after variable */
        const MemberRules& inherited_rules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inherited_rules.begin(), inherited_rules.end() );
        
        rules_set = true;
    }
    
    return move_member_rules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Move_MultiValueEventBirthDeath::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}



void Move_MultiValueEventBirthDeath::printValue(std::ostream &o) const
{
    
    o << "Move_MultiValueEventBirthDeath(";
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
void Move_MultiValueEventBirthDeath::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "x" )
    {
        x = var;
    }
    else if ( name == "ac" )
    {
        ac = var;
    }
//    else if ( name == "tune" )
//    {
//        tune = var;
//    }
    else
    {
        Move::setConstParameter(name, var);
    }
    
}
