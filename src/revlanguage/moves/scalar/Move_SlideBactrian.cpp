#include <cstddef>
#include <ostream>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RlBoolean.h"
#include "ContinuousStochasticNode.h"
#include "MetropolisHastingsMove.h"
#include "Move_SlideBactrian.h"
#include "Natural.h"
#include "Probability.h"
#include "Real.h"
#include "RealPos.h"
#include "RevObject.h"
#include "SlideBactrianProposal.h"
#include "SlideBactrianProposalContinuous.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "Move.h"
#include "RbBoolean.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlMove.h"
#include "StochasticNode.h"

namespace RevBayesCore { class Proposal; }


using namespace RevLanguage;

/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Move_SlideBactrian::Move_SlideBactrian() : Move()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the move.
 */
Move_SlideBactrian* Move_SlideBactrian::clone(void) const
{
    
    return new Move_SlideBactrian(*this);
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
void Move_SlideBactrian::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    RevBayesCore::Proposal *p = NULL;
    
    // now allocate a new sliding move
    double d = static_cast<const RealPos &>( sigma->getRevObject() ).getValue();
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    double r = static_cast<const RealPos &>( tune_target->getRevObject() ).getValue();
    size_t del = static_cast<const Natural &>( delay->getRevObject() ).getValue();
    RevBayesCore::TypedDagNode<double>* tmp = static_cast<const RealPos &>( x->getRevObject() ).getDagNode();
    RevBayesCore::ContinuousStochasticNode *n = dynamic_cast<RevBayesCore::ContinuousStochasticNode *>( tmp );
    if ( n != NULL )
    {
        p = new RevBayesCore::SlideBactrianProposalContinuous(n,d,r);
    }
    else
    {
        RevBayesCore::StochasticNode<double> *n2 = dynamic_cast<RevBayesCore::StochasticNode<double> *>( tmp );
        p = new RevBayesCore::SlideBactrianProposal(n2,d,r);
        
    }
    bool t = static_cast<const RlBoolean &>( tune->getRevObject() ).getValue();
    
    // finally create the internal move object
    //    value = new RevBayesCore::Move_SlideBactrian(n, d, t, w);
    
    value = new RevBayesCore::MetropolisHastingsMove(p,w,del,t);
    
}


/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Move_SlideBactrian::getClassType(void)
{
    
    static std::string rev_type = "Move_SlideBactrian";
    
    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Move_SlideBactrian::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_SlideBactrian::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "SlideBactrian";
    
    return c_name;
}


/**
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the SlideBactrian move are:
 * (1) the variable which must be a real.
 * (2) the tuning parameter lambda that defines the size of the proposal (positive real)
 * (3) a flag whether auto-tuning should be used.
 *
 * \return The member rules.
 */
const MemberRules& Move_SlideBactrian::getParameterRules(void) const
{
    
    static MemberRules slidingmove_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        slidingmove_member_rules.push_back( new ArgumentRule( "x"     , Real::getClassTypeSpec()     , "The variable on which this move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        slidingmove_member_rules.push_back( new ArgumentRule( "sigma" , RealPos::getClassTypeSpec()  , "The scale parameter.", ArgumentRule::BY_VALUE    , ArgumentRule::ANY       , new RealPos(1.0) ) );
        slidingmove_member_rules.push_back( new ArgumentRule( "tune"  , RlBoolean::getClassTypeSpec(), "Should we tune sigma during burnin?", ArgumentRule::BY_VALUE    , ArgumentRule::ANY       , new RlBoolean( true ) ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        slidingmove_member_rules.insert( slidingmove_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return slidingmove_member_rules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Move_SlideBactrian::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}



void Move_SlideBactrian::printValue(std::ostream &o) const
{
    
    o << "SlideBactrian(";
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
void Move_SlideBactrian::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "x" )
    {
        x = var;
    }
    else if ( name == "sigma" )
    {
        sigma = var;
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
