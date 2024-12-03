#include "AdaptiveReversibleJumpProposal.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MetropolisHastingsMove.h"
#include "Move_AdaptiveReversibleJumpSwitch.h"
#include "Natural.h"
#include "RbException.h"
#include "RealPos.h"
#include "TypeSpec.h"


using namespace RevLanguage;


Move_AdaptiveReversibleJumpSwitch::Move_AdaptiveReversibleJumpSwitch() : Move()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_AdaptiveReversibleJumpSwitch* Move_AdaptiveReversibleJumpSwitch::clone(void) const
{
    
    return new Move_AdaptiveReversibleJumpSwitch(*this);
}


void Move_AdaptiveReversibleJumpSwitch::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new vector-scale move
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    RevBayesCore::TypedDagNode<double>* tmp = static_cast<const Real &>( x->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<double> *sn = static_cast<RevBayesCore::StochasticNode<double> *>( tmp );
    long wbl = static_cast<const Natural &>( wait_before_learning->getRevObject() ).getValue();
    long wbu = static_cast<const Natural &>( wait_before_using->getRevObject() ).getValue();
    long ue  = static_cast<const Natural &>( update_every->getRevObject() ).getValue();
    RevBayesCore::AdaptiveReversibleJumpProposal *p = new RevBayesCore::AdaptiveReversibleJumpProposal(sn, wbl, wbu, ue);
    value = new RevBayesCore::MetropolisHastingsMove(p,w);

}


/** Get Rev type of object */
const std::string& Move_AdaptiveReversibleJumpSwitch::getClassType(void)
{
    
    static std::string rev_type = "Move_AdaptiveReversibleJumpSwitch";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& Move_AdaptiveReversibleJumpSwitch::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_AdaptiveReversibleJumpSwitch::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "AdaptiveRJSwitch";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Move_AdaptiveReversibleJumpSwitch::getParameterRules(void) const
{
    
    static MemberRules move_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        move_member_rules.push_back( new ArgumentRule( "x", Real::getClassTypeSpec(), "The variable on which this move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        move_member_rules.push_back( new ArgumentRule( "waitBeforeLearning", Natural::getClassTypeSpec(), "Number of tries before learning.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(1000) ) );
        move_member_rules.push_back( new ArgumentRule( "waitBeforeUsing", Natural::getClassTypeSpec(), "Number of tries before using.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(10000) ) );
        move_member_rules.push_back( new ArgumentRule( "updateEvery", Natural::getClassTypeSpec(), "How frequent to update.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(1) ) );

        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return move_member_rules;
}

/** Get type spec */
const TypeSpec& Move_AdaptiveReversibleJumpSwitch::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void Move_AdaptiveReversibleJumpSwitch::printValue(std::ostream &o) const
{
    
    o << "Move_AdaptiveReversibleJumpSwitch(";
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
void Move_AdaptiveReversibleJumpSwitch::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
    if ( name == "x" )
    {
        x = var;
    }
    else if ( name == "waitBeforeLearning" )
    {
        wait_before_learning = var;
    }
    else if ( name == "waitBeforeUsing" )
    {
        wait_before_using = var;
    }
    else if ( name == "updateEvery" )
    {
        update_every = var;
    }
    else
    {
        Move::setConstParameter(name, var);
    }
}


