#include <cstddef>
#include <ostream>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Move_EmpiricalTreeTopology.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlBoolean.h"
#include "RlTimeTree.h"
#include "RlTraceTree.h"
#include "TypeSpec.h"
#include "Move.h"
#include "RbBoolean.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlMove.h"
#include "StochasticNode.h"

namespace RevBayesCore { class TraceTree; }
namespace RevBayesCore { class Tree; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }


using namespace RevLanguage;

Move_EmpiricalTreeTopology::Move_EmpiricalTreeTopology() : Move()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_EmpiricalTreeTopology* Move_EmpiricalTreeTopology::clone(void) const
{
    
    return new Move_EmpiricalTreeTopology(*this);
}


void Move_EmpiricalTreeTopology::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new tree scale move
    RevBayesCore::TypedDagNode<RevBayesCore::Tree> *tmp = static_cast<const TimeTree &>( tree->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::Tree> *t = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree> *>( tmp );
    
    RevBayesCore::TraceTree* tt = &static_cast<const TraceTree &>( trace->getRevObject() ).getValue();

    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    double l = static_cast<const RealPos &>( alpha->getRevObject() ).getValue();
    bool tune = static_cast<const RlBoolean &>( tuning->getRevObject() ).getValue();
    
//    RevBayesCore::Proposal *p = new RevBayesCore::EmpiricalTreeTopologyProposal(t, tt, l);
//    value = new RevBayesCore::MetropolisHastingsMove(p, w, tune);
}


/** Get Rev type of object */
const std::string& Move_EmpiricalTreeTopology::getClassType(void)
{
    
    static std::string rev_type = "Move_EmpiricalTreeTopology";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& Move_EmpiricalTreeTopology::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_EmpiricalTreeTopology::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "EmpiricalTreeTopology";
    
    return c_name;
}


/** Return member rules */
const MemberRules& Move_EmpiricalTreeTopology::getParameterRules(void) const
{
    
    static MemberRules move_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        move_member_rules.push_back( new ArgumentRule( "tree"   , TimeTree::getClassTypeSpec() , "The tree variable the move operates on.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        move_member_rules.push_back( new ArgumentRule( "trace", TraceTree::getClassTypeSpec(), "The trace of tree samples.", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY ) );
        move_member_rules.push_back( new ArgumentRule( "alpha"  , RealPos::getClassTypeSpec()  , "The concentration factor of the proposal.", ArgumentRule::BY_VALUE    , ArgumentRule::ANY       , new RealPos( 1.0 ) ) );
        move_member_rules.push_back( new ArgumentRule( "tune"   , RlBoolean::getClassTypeSpec(), "Should we tune the scaling factor during burnin?", ArgumentRule::BY_VALUE    , ArgumentRule::ANY       , new RlBoolean( true ) ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return move_member_rules;
}

/** Get type spec */
const TypeSpec& Move_EmpiricalTreeTopology::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}



/** Get type spec */
void Move_EmpiricalTreeTopology::printValue(std::ostream &o) const
{
    
    o << "Move_EmpiricalTreeTopology(";
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
void Move_EmpiricalTreeTopology::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "tree" )
    {
        tree = var;
    }
    else if ( name == "trace" )
    {
        trace = var;
    }
    else if ( name == "alpha" )
    {
        alpha = var;
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
