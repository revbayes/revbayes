#include <cstddef>
#include <ostream>
#include <string>

#include <boost/optional.hpp>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "OptionRule.h"
#include "RlBoolean.h"
#include "ContinuousStochasticNode.h"
#include "Move_SliceSampling.h"
#include "Real.h"
#include "RealPos.h"
#include "RevObject.h"
#include "SliceSamplingMove.h"
#include "TypeSpec.h"
#include "Move.h"
#include "RbBoolean.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlMove.h"
#include "RlString.h"

using boost::optional;

namespace RevBayesCore { template <class valueType> class TypedDagNode; }


using namespace RevLanguage;

/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Move_SliceSampling::Move_SliceSampling() : Move()
{

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the move.
 */
Move_SliceSampling* Move_SliceSampling::clone(void) const
{

	return new Move_SliceSampling(*this);
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
void Move_SliceSampling::constructInternalObject( void )
{
    // we free the memory first
    delete value;

    // now allocate a new sliding move
    double window_ = static_cast<const RealPos &>( window->getRevObject() ).getValue();
    double weight_ = static_cast<const RealPos &>( weight->getRevObject() ).getValue();

    boost::optional<double> lower_bound_;
    boost::optional<double> upper_bound_;
    // If the x variable is a RealPos, the lower bound needs to be at least zero.
    if (dynamic_cast<const RealPos *>( & x->getRevObject() ))
    {
        if (lower_bound_)
            lower_bound_ = std::max(0.0, *lower_bound_);
        else
            lower_bound_ = 0.0;
    }

    RevBayesCore::TypedDagNode<double>* tmp = static_cast<const Real &>( x->getRevObject() ).getDagNode();
    RevBayesCore::ContinuousStochasticNode *node_ = static_cast<RevBayesCore::ContinuousStochasticNode *>( tmp );
    bool tune_ = static_cast<const RlBoolean &>( tune->getRevObject() ).getValue();

    auto search_method_ = dynamic_cast<const RlString&>( search_method->getRevObject() ).getValue();
    typedef RevBayesCore::SliceSamplingMove::BoundarySearchMethod search_method_t;
    search_method_t search_method__;
    if (search_method_ == "doubling")
        search_method__ = RevBayesCore::SliceSamplingMove::BoundarySearchMethod::search_doubling;
    else if (search_method_ == "stepping_out")
        search_method__ = RevBayesCore::SliceSamplingMove::BoundarySearchMethod::search_stepping_out;
    else
        throw RbException("mvSlice: search_method should be \"doubling\" or \"stepping out\"");

    // finally create the internal move object

    value = new RevBayesCore::SliceSamplingMove(node_ , lower_bound_, upper_bound_, window_, weight_ , search_method__, tune_);
}


/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Move_SliceSampling::getClassType(void)
{

    static std::string rev_type = "Move_SliceScampling";

	return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Move_SliceSampling::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );

	return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_SliceSampling::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "Slice";

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
const MemberRules& Move_SliceSampling::getParameterRules(void) const
{

    static MemberRules move_member_rules;
    static bool rules_set = false;

    if ( !rules_set )
    {
        move_member_rules.push_back( new ArgumentRule( "x"     ,  Real::getClassTypeSpec()    , "The variable on which this move operates", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        move_member_rules.push_back( new ArgumentRule( "window", RealPos::getClassTypeSpec()  , "The window (steps-size) of proposals.", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new RealPos(1.0) ) );
        move_member_rules.push_back( new ArgumentRule( "tune"  , RlBoolean::getClassTypeSpec(), "Should we tune the move during burnin?", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new RlBoolean( true ) ) );

        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );

        std::vector<std::string> optionsMethod;
        optionsMethod.push_back( "stepping_out" );
        optionsMethod.push_back( "doubling" );
        move_member_rules.push_back( new OptionRule( "search_method", new RlString("doubling"), optionsMethod, "The method used to find the slice boundaries." ) );

        rules_set = true;
    }

    return move_member_rules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Move_SliceSampling::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}



void Move_SliceSampling::printValue(std::ostream &o) const {

    o << "SliceSampling(";
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
void Move_SliceSampling::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{

    if ( name == "x" )
    {
        x = var;
    }
    else if ( name == "window" )
    {
        window = var;
    }
    else if ( name == "tune" )
    {
        tune = var;
    }
    else if ( name == "search_method" )
    {
        search_method = var;
    }
    else
    {
        Move::setConstParameter(name, var);
    }

}
