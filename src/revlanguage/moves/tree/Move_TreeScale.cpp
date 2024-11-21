#include <cstddef>
#include <ostream>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MetropolisHastingsMove.h"
#include "ModelVector.h"
#include "Move_TreeScale.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlBoolean.h"
#include "RlTimeTree.h"
#include "TreeScaleProposal.h"
#include "TypeSpec.h"
#include "Move.h"
#include "RbBoolean.h"
#include "RevNullObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlMove.h"
#include "StochasticNode.h"


using namespace RevLanguage;

Move_TreeScale::Move_TreeScale() : Move()
{

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_TreeScale* Move_TreeScale::clone(void) const
{

	return new Move_TreeScale(*this);
}


void Move_TreeScale::constructInternalObject( void )
{
    // we free the memory first
    delete value;

    // now allocate a new tree scale move

    // get the tree(s) variable
    // we either expect to receive a stochastic variable on a single tree (e.g., the species tree), or on a vector of trees (e.g., gene trees)
    // to avoid re-implementing the move for a single scalar and vector, we allow for both single values and vectors
    // however, we need to check to make sure that our type conversions stay proper
    RevBayesCore::StochasticNode<RevBayesCore::Tree> *t = NULL;
    RevBayesCore::StochasticNode< RevBayesCore::RbVector<RevBayesCore::Tree> > *vec_t = NULL;
    if ( tree->getRevObject().isType( TimeTree::getClassTypeSpec() ) )
    {
        RevBayesCore::TypedDagNode<RevBayesCore::Tree> *tmp = static_cast<const TimeTree &>( tree->getRevObject() ).getDagNode();
        t = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree> *>( tmp );
    }
    else if ( tree->getRevObject().isType( ModelVector<TimeTree>::getClassTypeSpec() ) )
    {
        RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::Tree> > *tmp = static_cast<const ModelVector<TimeTree> &>( tree->getRevObject() ).getDagNode();
        vec_t = static_cast<RevBayesCore::StochasticNode< RevBayesCore::RbVector<RevBayesCore::Tree> > *>( tmp );
    }
    else
    {
        throw RbException() << "Wrong tree type '" <<  tree->getRevObject().getType() << "'.";
    }

    RevBayesCore::StochasticNode<double> *ra = NULL;
    if ( rootAge != NULL && rootAge->getRevObject() != RevNullObject::getInstance() )
    {
        RevBayesCore::TypedDagNode<double> *tmp = static_cast<const RealPos &>( rootAge->getRevObject() ).getDagNode();
        ra = static_cast<RevBayesCore::StochasticNode<double> *>( tmp );
    }
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    double l = static_cast<const RealPos &>( delta->getRevObject() ).getValue();
    bool tune = static_cast<const RlBoolean &>( tuning->getRevObject() ).getValue();

    RevBayesCore::Proposal *p = new RevBayesCore::TreeScaleProposal(t, vec_t, ra, l);
    value = new RevBayesCore::MetropolisHastingsMove(p, w, tune);
}


/** Get Rev type of object */
const std::string& Move_TreeScale::getClassType(void)
{

    static std::string rev_type = "Move_TreeScale";

	return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& Move_TreeScale::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );

	return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_TreeScale::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "TreeScale";

    return c_name;
}


/** Return member rules */
const MemberRules& Move_TreeScale::getParameterRules(void) const
{

    static MemberRules move_member_rules;
    static bool rules_set = false;

    if ( !rules_set )
    {
        std::vector<TypeSpec> tree_var_types;
        tree_var_types.push_back( TimeTree::getClassTypeSpec() );
        tree_var_types.push_back( ModelVector<TimeTree>::getClassTypeSpec() );
        move_member_rules.push_back( new ArgumentRule( "tree"   , tree_var_types               , "The tree on which this move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        move_member_rules.push_back( new ArgumentRule( "rootAge", RealPos::getClassTypeSpec()  , "The root age variable.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC, NULL ) );
        move_member_rules.push_back( new ArgumentRule( "delta"  , RealPos::getClassTypeSpec()  , "The scaling factor (strength) of the proposal.", ArgumentRule::BY_VALUE    , ArgumentRule::ANY       , new RealPos( 1.0 ) ) );
        move_member_rules.push_back( new ArgumentRule( "tune"   , RlBoolean::getClassTypeSpec(), "Should we tune the scaling factor during burnin?", ArgumentRule::BY_VALUE    , ArgumentRule::ANY       , new RlBoolean( true ) ) );

        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );

        rules_set = true;
    }

    return move_member_rules;
}

/** Get type spec */
const TypeSpec& Move_TreeScale::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}



/** Get type spec */
void Move_TreeScale::printValue(std::ostream &o) const
{

    o << "Move_TreeScale(";
    if (tree != NULL)
    {
        o << tree->getName();
    }
    else
    {
        o << "?";
    }
    o << ")";
}


/** Set a NearestNeighborInterchange variable */
void Move_TreeScale::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{

    if ( name == "tree" )
    {
        tree = var;
    }
    else if ( name == "rootAge" )
    {
        rootAge = var;
    }
    else if ( name == "delta" )
    {
        delta = var;
    }
    else if ( name == "tune" )
    {
        tuning = var;
    }
    else
    {
        Move::setConstParameter(name, var);
    }

}
