#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RlBoolean.h"
#include "EllipticalSliceSamplingLognormalIIDMove.h"
#include "ModelVector.h"
#include "Move_EllipticalSliceSamplingLognormalIID.h"
#include "RbException.h"
#include "RbConstants.h"
#include "Real.h"
#include "RealPos.h"
#include "RevObject.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "DagNode.h"
#include "ModelObject.h"
#include "Move.h"
#include "RbBoolean.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlMove.h"
#include "StochasticNode.h"

namespace RevBayesCore { template <class valueType> class RbVector; }


using namespace RevLanguage;

/**
 * Default constructor.
 * 
 * The default constructor does nothing except allocating the object.
 */
Move_EllipticalSliceSamplingLognormalIID::Move_EllipticalSliceSamplingLognormalIID() : Move() 
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the move. 
 */
Move_EllipticalSliceSamplingLognormalIID* Move_EllipticalSliceSamplingLognormalIID::clone(void) const 
{
    
	return new Move_EllipticalSliceSamplingLognormalIID(*this);
}


/**
 * Create a new internal move object.
 *
 * This function simply dynamically allocates a new internal move object that is 
 * associated with the variable (DAG-node). The internal move object is created by calling its
 * constructor and passing the move-parameters (the variable and other parameters) as arguments of the 
 * constructor. The move constructor takes care of the proper hook-ups.
 *
 * \return A new internal distribution object.
 */
void Move_EllipticalSliceSamplingLognormalIID::constructInternalObject( void ) 
{
    // we free the memory first
    delete value;
    
    // now allocate a new sliding move
    double window_ = static_cast<const RealPos &>( window->getRevObject() ).getValue();
    double weight_ = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    
    // Add the field parameters
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* tmp = static_cast<const ModelVector<RealPos> &>( x->getRevObject() ).getDagNode();
    std::vector<const RevBayesCore::DagNode*> p = tmp->getParents();
    std::vector< RevBayesCore::StochasticNode<double> *> n;
    for (std::vector<const RevBayesCore::DagNode*>::const_iterator it = p.begin(); it != p.end(); ++it)
    {
        const RevBayesCore::StochasticNode<double> *the_node = dynamic_cast< const RevBayesCore::StochasticNode<double>* >( *it );
        if ( the_node != NULL )
        {
            n.push_back( const_cast< RevBayesCore::StochasticNode<double>* >( the_node ) );
        }
        else
        {
            throw RbException("Could not create a mvEllipticalSliceSamplingLognormalIID because the node isn't a vector of stochastic nodes.");
        }
    }
    
    RevBayesCore::TypedDagNode< double >* m = static_cast<const Real &>( mu->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< double >* s = static_cast<const RealPos &>( sigma->getRevObject() ).getDagNode();

    bool tune_ = static_cast<const RlBoolean &>( tune->getRevObject() ).getValue();
    
    // finally create the internal move object
    
    value = new RevBayesCore::EllipticalSliceSamplingLognormalIIDMove(n, m, s, window_, weight_ ,tune_);
}


/**
 * Get Rev type of object 
 *
 * \return The class' name.
 */
const std::string& Move_EllipticalSliceSamplingLognormalIID::getClassType(void) 
{ 
    
    static std::string rev_type = "Move_EllipticalSliceSamplingLognormalIID";
    
	return rev_type; 
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Move_EllipticalSliceSamplingLognormalIID::getClassTypeSpec(void) 
{ 
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_EllipticalSliceSamplingLognormalIID::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "EllipticalSliceSamplingLognormalIID";
    
    return c_name;
}


/** 
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the scale move are:
 * (1) the variable which must be a positive real.
 * (2) the tuning parameter window that defines the size of the proposal (positive real)
 * (3) a flag whether auto-tuning should be used. 
 *
 * \return The member rules.
 */
const MemberRules& Move_EllipticalSliceSamplingLognormalIID::getParameterRules(void) const 
{
    
    static MemberRules move_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set ) 
    {
        move_member_rules.push_back( new ArgumentRule( "x"     , ModelVector<RealPos>::getClassTypeSpec(), "The variable on which this move operates", ArgumentRule::BY_REFERENCE, ArgumentRule::DETERMINISTIC ) );
        move_member_rules.push_back( new ArgumentRule( "mu",     Real::getClassTypeSpec(),                 "Log-mean parameter of distribution.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        move_member_rules.push_back( new ArgumentRule( "sigma" , RealPos::getClassTypeSpec(),              "Log-sd parameter of distribution.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        move_member_rules.push_back( new ArgumentRule( "window", RealPos::getClassTypeSpec(),              "The anglular window (steps-size) of proposals, a value in the interval (0,2*pi].", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new RealPos(RbConstants::TwoPI) ) );
        move_member_rules.push_back( new ArgumentRule( "tune"  , RlBoolean::getClassTypeSpec(),            "Should we tune the move during burnin?", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new RlBoolean( true ) ) );

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
const TypeSpec& Move_EllipticalSliceSamplingLognormalIID::getTypeSpec( void ) const 
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}



void Move_EllipticalSliceSamplingLognormalIID::printValue(std::ostream &o) const {
    
    o << "Move_EllipticalSliceSamplingLognormalIID(";
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
void Move_EllipticalSliceSamplingLognormalIID::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) 
{
    
    if ( name == "x" ) 
    {
        x = var;
    }
    else if ( name == "mu" )
    {
        mu = var;
    }
    else if ( name == "sigma" )
    {
        sigma = var;
    }
    else if ( name == "window" )
    {
        window = var;
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
