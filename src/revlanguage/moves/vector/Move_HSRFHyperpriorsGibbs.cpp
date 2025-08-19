#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "HSRFHyperpriorsGibbsMove.h"
#include "HSRFOrder2HyperpriorsGibbsMove.h"
#include "ModelVector.h"
#include "Move_HSRFHyperpriorsGibbs.h"
#include "Natural.h"
#include "Probability.h"
#include "RbException.h"
#include "Real.h"
#include "RealPos.h"
#include "RevObject.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "DagNode.h"
#include "ModelObject.h"
#include "Move.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlMove.h"
#include "StochasticNode.h"

namespace RevBayesCore { template <class valueType> class RbVector; }


using namespace RevLanguage;

Move_HSRFHyperpriorsGibbs::Move_HSRFHyperpriorsGibbs() : Move()
{

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_HSRFHyperpriorsGibbs* Move_HSRFHyperpriorsGibbs::clone(void) const
{

    return new Move_HSRFHyperpriorsGibbs(*this);
}


void Move_HSRFHyperpriorsGibbs::constructInternalObject( void )
{
    // we free the memory first
    delete value;

    // now allocate a new move
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();

    double p = static_cast<const RealPos &>( prop_global->getRevObject() ).getValue();

    RevBayesCore::TypedDagNode<double>* tmp = static_cast<const RealPos &>( gs->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<double> *g = static_cast< RevBayesCore::StochasticNode<double>* >( tmp );

    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* tmp1 = static_cast<const ModelVector<RealPos> &>( ls->getRevObject() ).getDagNode();
    std::vector<const RevBayesCore::DagNode*> p1 = tmp1->getParents();
    std::vector< RevBayesCore::StochasticNode<double> *> l;
    for (std::vector<const RevBayesCore::DagNode*>::const_iterator it = p1.begin(); it != p1.end(); ++it)
    {
        const RevBayesCore::StochasticNode<double> *the_node = dynamic_cast< const RevBayesCore::StochasticNode<double>* >( *it );
        if ( the_node != NULL )
        {
            l.push_back( const_cast< RevBayesCore::StochasticNode<double>* >( the_node ) );
        }
        else
        {
            throw RbException("Could not create a mvHSRFHyperpriorsGibbs because the node isn't a vector of stochastic nodes.");
        }
    }

    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* tmp2 = static_cast<const ModelVector<RealPos> &>( normals->getRevObject() ).getDagNode();
    std::vector<const RevBayesCore::DagNode*> p2 = tmp2->getParents();
    std::vector< RevBayesCore::StochasticNode<double> *> n;
    for (std::vector<const RevBayesCore::DagNode*>::const_iterator it = p2.begin(); it != p2.end(); ++it)
    {
        const RevBayesCore::StochasticNode<double> *the_node = dynamic_cast< const RevBayesCore::StochasticNode<double>* >( *it );
        if ( the_node != NULL )
        {
            n.push_back( const_cast< RevBayesCore::StochasticNode<double>* >( the_node ) );
        }
        else
        {
            throw RbException("Could not create a mvHSRFHyperpriorGibbs because \"normals\" isn't a vector of stochastic nodes.");
        }
    }

    double z = static_cast<const RealPos &>( zeta->getRevObject() ).getValue();

    size_t o = static_cast<const Natural &>( order->getRevObject() ).getValue();
    if (o == 1)
    {
      value = new RevBayesCore::HSRFHyperpriorsGibbsMove(g, l, n, p, z, w);
    }
    else if (o == 2)
    {
      value = new RevBayesCore::HSRFOrder2HyperpriorsGibbsMove(g, l, n, p, z, w);
    }
    else
    {
      throw(RbException("Only valid options for option \"order\" in mvHSRFHyperpriorGibbs are 1 and 2."));
    }

}


/** Get Rev type of object */
const std::string& Move_HSRFHyperpriorsGibbs::getClassType(void)
{

    static std::string rev_type = "Move_HSRFHyperpriorsGibbs";

    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& Move_HSRFHyperpriorsGibbs::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_HSRFHyperpriorsGibbs::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "HSRFHyperpriorsGibbs";

    return c_name;
}



/** Return member rules (no members) */
const MemberRules& Move_HSRFHyperpriorsGibbs::getParameterRules(void) const
{

    static MemberRules move_member_rules;
    static bool rules_set = false;

    if ( !rules_set )
    {

        move_member_rules.push_back( new ArgumentRule( "gs"            , RealPos::getClassTypeSpec(),              "The global scale variable on which this move operates", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        move_member_rules.push_back( new ArgumentRule( "ls"            , ModelVector<RealPos>::getClassTypeSpec(), "The local scales variables on which this move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::DETERMINISTIC ) );
        move_member_rules.push_back( new ArgumentRule( "normals"       , ModelVector<Real>::getClassTypeSpec(),    "The vector of Normal RVs defining the field", ArgumentRule::BY_REFERENCE, ArgumentRule::DETERMINISTIC ) );
        move_member_rules.push_back( new ArgumentRule( "zeta"          , RealPos::getClassTypeSpec()  ,            "The value controlling the shrinkage of the field.", ArgumentRule::BY_VALUE    , ArgumentRule::ANY ) );
        move_member_rules.push_back( new ArgumentRule( "propGlobalOnly", Probability::getClassTypeSpec()  ,        "The proportion of all of these moves that only move the global scale.", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new Probability(0.0) ) );
        move_member_rules.push_back( new ArgumentRule( "order"         , Natural::getClassTypeSpec()  ,            "The order of this HSRF model, first (1) or second (2). Defaults to first order.", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new Natural(1) ) );

        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );

        rules_set = true;
    }

    return move_member_rules;
}

/** Get type spec */
const TypeSpec& Move_HSRFHyperpriorsGibbs::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}


/** Get type spec */
void Move_HSRFHyperpriorsGibbs::printValue(std::ostream &o) const
{

    o << "Move_HSRFHyperpriorsGibbsMove(";
    if (gs != NULL && ls != NULL)
    {
        o << gs->getName();
        o << ",";
        o << ls->getName();
    }
    else
    {
        o << "?";
    }
    o << ")";
}


/** Set a member variable */
void Move_HSRFHyperpriorsGibbs::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {

    if ( name == "gs" )
    {
        gs = var;
    }
    else if ( name == "ls" )
    {
        ls = var;
    }

    else if ( name == "normals" )
    {
        normals = var;
    }
    else if ( name == "zeta" )
    {
        zeta = var;
    }
    else if ( name == "propGlobalOnly" )
    {
        prop_global = var;
    }
    else if ( name == "order" )
    {
        order = var;
    }
    else
    {
        Move::setConstParameter(name, var);
    }

}
