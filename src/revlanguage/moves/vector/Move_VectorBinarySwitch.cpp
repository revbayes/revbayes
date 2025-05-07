#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MetropolisHastingsMove.h"
#include "ModelVector.h"
#include "Move_VectorBinarySwitch.h"
#include "Natural.h"
#include "RbException.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlBoolean.h"
#include "Probability.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "VectorBinarySwitchProposal.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "Move.h"
#include "RbBoolean.h"
#include "RbVector.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "RlMove.h"
#include "StochasticNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

namespace RevBayesCore { class Proposal; }


using namespace RevLanguage;

Move_VectorBinarySwitch::Move_VectorBinarySwitch() : Move()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_VectorBinarySwitch* Move_VectorBinarySwitch::clone(void) const
{
    
    return new Move_VectorBinarySwitch(*this);
}


void Move_VectorBinarySwitch::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new vector-slide move
    const RevBayesCore::RbVector<long> &e = static_cast<const ModelVector<Natural> &>( indices->getRevObject() ).getValue();
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    double p = static_cast<const Probability &>( switch_probability->getRevObject() ).getValue();
    
    
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<long> >* tmp = static_cast<const ModelVector<Natural> &>( x->getRevObject() ).getDagNode();
//    RevBayesCore::StochasticNode<RevBayesCore::RbVector<long> > *n = static_cast<RevBayesCore::StochasticNode<RevBayesCore::RbVector<long> > *>( tmp );
    
    std::vector<const RevBayesCore::DagNode*> dag_nodes = tmp->getParents();
    std::vector< RevBayesCore::StochasticNode<long>* > stoch_nodes;
    for (std::vector<const RevBayesCore::DagNode*>::const_iterator it = dag_nodes.begin(); it != dag_nodes.end(); ++it)
    {
        const RevBayesCore::StochasticNode<long> *the_node = dynamic_cast< const RevBayesCore::StochasticNode<long>* >( *it );
        if ( the_node != NULL )
        {
            stoch_nodes.push_back( const_cast< RevBayesCore::StochasticNode<long>* >( the_node ) );
        }
        else
        {
            throw RbException("Could not create a mvShrinkExpand because the node isn't a vector of stochastic nodes.");
        }
    }

    
    
    
    bool t = static_cast<const RlBoolean &>( tune->getRevObject() ).getValue();
    
    RevBayesCore::Proposal *prop = new RevBayesCore::VectorBinarySwitchProposal(stoch_nodes,e,p);
    value = new RevBayesCore::MetropolisHastingsMove(prop,w,t);
}


/** Get Rev type of object */
const std::string& Move_VectorBinarySwitch::getClassType(void)
{
    
    static std::string rev_type = "Move_VectorBinarySwitch";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& Move_VectorBinarySwitch::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_VectorBinarySwitch::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "VectorBinarySwitch";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Move_VectorBinarySwitch::getParameterRules(void) const
{
    
    static MemberRules move_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        move_member_rules.push_back( new ArgumentRule( "x"        , ModelVector<Natural>::getClassTypeSpec()   , "The variable on which this move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::DETERMINISTIC ) );
        move_member_rules.push_back( new ArgumentRule( "p"        , Probability::getClassTypeSpec()         , "The probability of applying the move to each element.", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY, new Probability(1.0) ) );
        move_member_rules.push_back( new ArgumentRule( "elements" , ModelVector<Natural>::getClassTypeSpec(), "The indices of elements. All are used if this is empty.", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY, new ModelVector<Natural>() ) );
        move_member_rules.push_back( new ArgumentRule( "tune"     , RlBoolean::getClassTypeSpec()           , "Should we tune the window size during burnin?", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new RlBoolean( true ) ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return move_member_rules;
}

/** Get type spec */
const TypeSpec& Move_VectorBinarySwitch::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void Move_VectorBinarySwitch::printValue(std::ostream &o) const
{
    
    o << "Move_VectorBinarySwitch(";
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
void Move_VectorBinarySwitch::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
    if ( name == "x" )
    {
        x = var;
    }
    else if (name == "p")
    {
        switch_probability = var;
    }
    else if ( name == "elements" )
    {
        indices = var;
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
