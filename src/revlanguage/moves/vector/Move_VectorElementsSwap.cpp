#include <stddef.h>
#include <ostream>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RlBoolean.h"
#include "MetropolisHastingsMove.h"
#include "ModelVector.h"
#include "Move_VectorElementsSwap.h"
#include "Natural.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlSimplex.h"
#include "VectorElementsSwapProposal.h"
#include "TypeSpec.h"
#include "Move.h"
#include "RbBoolean.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlMove.h"
#include "StochasticNode.h"

namespace RevBayesCore { class Proposal; }
namespace RevBayesCore { template <class valueType> class RbVector; }


using namespace RevLanguage;

Move_VectorElementsSwap::Move_VectorElementsSwap() : Move()
{
    // first, the argument rules
    ArgumentRules* addIndexRules = new ArgumentRules();
    
    // next, set the specific arguments
    addIndexRules->push_back( new ArgumentRule( "elementIndex", Natural::getClassTypeSpec(), "The index of the element to add to the swapping set.", ArgumentRule::BY_VALUE , ArgumentRule::ANY, new Natural( 1 ) ) );
    
    // finally, create the methods
    methods.addFunction( new MemberProcedure( "addIndex", RlUtils::Void, addIndexRules) );
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_VectorElementsSwap* Move_VectorElementsSwap::clone(void) const
{
    
    return new Move_VectorElementsSwap(*this);
}


void Move_VectorElementsSwap::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new move
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
            throw RbException("Could not create a mvElementScale because the node isn't a vector of stochastic nodes.");
        }
    }
    
    RevBayesCore::Proposal *prop = new RevBayesCore::VectorElementsSwapProposal(n);
    value = new RevBayesCore::MetropolisHastingsMove(prop, w);
    
}


RevPtr<RevVariable> Move_VectorElementsSwap::executeMethod(const std::string& name, const std::vector<Argument>& args, bool &found)
{
    
    if ( name == "addIndex" )
    {
        found = true;
        
        size_t idx = static_cast<const Natural &>( args[0].getVariable()->getRevObject() ).getValue();
        
        RevBayesCore::MetropolisHastingsMove *m = static_cast<RevBayesCore::MetropolisHastingsMove*>(this->value);
        RevBayesCore::VectorElementsSwapProposal &prop = static_cast<RevBayesCore::VectorElementsSwapProposal&>( m->getProposal() );
        
        RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* tmp = static_cast<const ModelVector<RealPos> &>( x->getRevObject() ).getDagNode();
        
        if ( idx >= 1 && idx <= tmp->getValue().size())
        {
            prop.addIndex( idx - 1 );
        }
        else
        {
            std::stringstream ss_err;
            ss_err << "Index out of bounds: The vector of size " << tmp->getValue().size() << " does not have an element for index " << idx << ".";
            throw RbException(ss_err.str());
        }
        
        return NULL;
    }
    
    return Move::executeMethod( name, args, found );
}


/** Get Rev type of object */
const std::string& Move_VectorElementsSwap::getClassType(void)
{
    
    static std::string rev_type = "Move_VectorElementsSwap";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Move_VectorElementsSwap::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_VectorElementsSwap::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "VectorElementsSwap";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Move_VectorElementsSwap::getParameterRules(void) const
{
    
    static MemberRules move_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        move_member_rules.push_back( new ArgumentRule( "x"    , ModelVector<RealPos>::getClassTypeSpec(), "The variable on which the move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::DETERMINISTIC ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return move_member_rules;
}

/** Get type spec */
const TypeSpec& Move_VectorElementsSwap::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void Move_VectorElementsSwap::printValue(std::ostream &o) const
{
    
    o << "Move_VectorElementsSwap(";
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
void Move_VectorElementsSwap::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "x" )
    {
        x = var;
    }
    else
    {
        Move::setConstParameter(name, var);
    }
}
