#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MetropolisHastingsMove.h"
#include "ModelVector.h"
#include "Move_ShrinkExpandScale.h"
#include "RbException.h"
#include "Real.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlBoolean.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "ShrinkExpandScaleProposal.h"
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

Move_ShrinkExpandScale::Move_ShrinkExpandScale() : Move()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_ShrinkExpandScale* Move_ShrinkExpandScale::clone(void) const
{
    
    return new Move_ShrinkExpandScale(*this);
}


void Move_ShrinkExpandScale::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new vector-scale move
    double l = static_cast<const RealPos &>( lambda->getRevObject() ).getValue();
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    
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
            throw RbException("Could not create a mvShrinkExpandScale because the node isn't a vector of stochastic nodes.");
        }
    }
    
    // get the standard deviation if we want to scale that one too
    RevBayesCore::StochasticNode<double> *s = NULL;
    if ( sd != NULL && sd->getRevObject() != RevNullObject::getInstance() )
    {
        RevBayesCore::TypedDagNode<double> *tmp = static_cast<const RealPos &>( sd->getRevObject() ).getDagNode();
        s = static_cast<RevBayesCore::StochasticNode<double> *>( tmp );
    }
    
    // get the tuning
    bool t = static_cast<const RlBoolean &>( tune->getRevObject() ).getValue();
    
    RevBayesCore::Proposal *prop = new RevBayesCore::ShrinkExpandScaleProposal(n,s,l);
    value = new RevBayesCore::MetropolisHastingsMove(prop,w,t);
    
}


/** Get Rev type of object */
const std::string& Move_ShrinkExpandScale::getClassType(void)
{
    
    static std::string rev_type = "Move_ShrinkExpandScale";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& Move_ShrinkExpandScale::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_ShrinkExpandScale::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "ShrinkExpandScale";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Move_ShrinkExpandScale::getParameterRules(void) const
{
    
    static MemberRules move_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        move_member_rules.push_back( new ArgumentRule( "x"     , ModelVector<RealPos>::getClassTypeSpec(), "The variable on which the move operates."        , ArgumentRule::BY_REFERENCE, ArgumentRule::DETERMINISTIC ) );
        move_member_rules.push_back( new ArgumentRule( "sd"    , RealPos::getClassTypeSpec()             , "The standard deviation parameter if available."  , ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC, NULL ) );
        move_member_rules.push_back( new ArgumentRule( "lambda", RealPos::getClassTypeSpec()             , "The scaling factor (strength) of the proposal."  , ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new Real(1.0) ) );
        move_member_rules.push_back( new ArgumentRule( "tune"  , RlBoolean::getClassTypeSpec()           , "Should we tune the scaling factor during burnin?", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new RlBoolean( true ) ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inherited_rules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inherited_rules.begin(), inherited_rules.end() );
        
        rules_set = true;
    }
    
    return move_member_rules;
}

/** Get type spec */
const TypeSpec& Move_ShrinkExpandScale::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void Move_ShrinkExpandScale::printValue(std::ostream &o) const
{
    
    o << "Move_ShrinkExpandScale(";
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
void Move_ShrinkExpandScale::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "x" )
    {
        x = var;
    }
    else if ( name == "sd" )
    {
        sd = var;
    }
    else if ( name == "lambda" )
    {
        lambda = var;
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
