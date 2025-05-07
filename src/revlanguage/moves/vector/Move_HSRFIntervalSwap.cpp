#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MetropolisHastingsMove.h"
#include "ModelVector.h"
#include "Move_HSRFIntervalSwap.h"
#include "RbException.h"
#include "Real.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlBoolean.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "HSRFIntervalSwapProposal.h"
#include "DagNode.h"
#include "ModelObject.h"
#include "Move.h"
#include "RbBoolean.h"
#include "RevNullObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlMove.h"
#include "StochasticNode.h"

namespace RevBayesCore { class Proposal; }
namespace RevBayesCore { template <class valueType> class RbVector; }


using namespace RevLanguage;

Move_HSRFIntervalSwap::Move_HSRFIntervalSwap() : Move()
{

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_HSRFIntervalSwap* Move_HSRFIntervalSwap::clone(void) const
{

    return new Move_HSRFIntervalSwap(*this);
}


void Move_HSRFIntervalSwap::constructInternalObject( void )
{
    // we free the memory first
    delete value;

    // now allocate a new vector-scale move
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();

    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* tmp1 = static_cast<const ModelVector<RealPos> &>( sigmas->getRevObject() ).getDagNode();
    std::vector<const RevBayesCore::DagNode*> p1 = tmp1->getParents();
    std::vector< RevBayesCore::StochasticNode<double> *> s;
    for (std::vector<const RevBayesCore::DagNode*>::const_iterator it = p1.begin(); it != p1.end(); ++it)
    {
        const RevBayesCore::StochasticNode<double> *the_node = dynamic_cast< const RevBayesCore::StochasticNode<double>* >( *it );
        if ( the_node != NULL )
        {
            s.push_back( const_cast< RevBayesCore::StochasticNode<double>* >( the_node ) );
        }
        else
        {
            throw RbException("Could not create a mvHSRFIntervalSwap because \"sigmas\" isn't a vector of stochastic nodes.");
        }
    }

    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* tmp2 = static_cast<const ModelVector<Real> &>( deltas->getRevObject() ).getDagNode();
    std::vector<const RevBayesCore::DagNode*> p2 = tmp2->getParents();
    std::vector< RevBayesCore::StochasticNode<double> *> d;
    for (std::vector<const RevBayesCore::DagNode*>::const_iterator it = p2.begin(); it != p2.end(); ++it)
    {
        const RevBayesCore::StochasticNode<double> *the_node = dynamic_cast< const RevBayesCore::StochasticNode<double>* >( *it );
        if ( the_node != NULL )
        {
            d.push_back( const_cast< RevBayesCore::StochasticNode<double>* >( the_node ) );
        }
        else
        {
            throw RbException("Could not create a mvHSRFIntervalSwap because \"deltas\" isn't a vector of stochastic nodes.");
        }
    }

    // get the tuning
    bool t = static_cast<const RlBoolean &>( tune->getRevObject() ).getValue();

    RevBayesCore::Proposal *prop = new RevBayesCore::HSRFIntervalSwapProposal(d,s);
    value = new RevBayesCore::MetropolisHastingsMove(prop,w,t);

}


/** Get Rev type of object */
const std::string& Move_HSRFIntervalSwap::getClassType(void)
{

    static std::string rev_type = "Move_HSRFIntervalSwap";

    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& Move_HSRFIntervalSwap::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_HSRFIntervalSwap::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "HSRFIntervalSwap";

    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Move_HSRFIntervalSwap::getParameterRules(void) const
{

    static MemberRules move_member_rules;
    static bool rules_set = false;

    if ( !rules_set )
    {

        move_member_rules.push_back( new ArgumentRule( "deltas", ModelVector<Real>::getClassTypeSpec(), "The HSMRF delta parameters on which the move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::DETERMINISTIC ) );
        move_member_rules.push_back( new ArgumentRule( "sigmas", ModelVector<RealPos>::getClassTypeSpec(), "The HSMRF sigma parameters on which the move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::DETERMINISTIC ) );
        move_member_rules.push_back( new ArgumentRule( "tune"  , RlBoolean::getClassTypeSpec()        , "Should we tune the scaling factor during burnin?", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new RlBoolean( true ) ) );

        /* Inherit weight from Move, put it after variable */
        const MemberRules& inherited_rules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inherited_rules.begin(), inherited_rules.end() );

        rules_set = true;
    }

    return move_member_rules;
}

/** Get type spec */
const TypeSpec& Move_HSRFIntervalSwap::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}


/** Get type spec */
void Move_HSRFIntervalSwap::printValue(std::ostream &o) const
{

    o << "Move_HSRFIntervalSwap(";
    if (deltas != NULL)
    {
        o << deltas->getName();
    }
    else
    {
        o << "?";
    }
    o << ",";
    if (sigmas != NULL)
    {
        o << sigmas->getName();
    }
    else
    {
        o << "?";
    }
    o << ")";

}


/** Set a member variable */
void Move_HSRFIntervalSwap::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{

    if ( name == "deltas" )
    {
        deltas = var;
    }
    else if ( name == "sigmas" )
    {
        sigmas = var;
    }
    else if ( name == "tune" )
    {
        tune = var;
    }
    else {
        Move::setConstParameter(name, var);
    }

}
