#include <cstddef>
#include <ostream>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MetropolisHastingsMove.h"
#include "ModelVector.h"
#include "Move_NodeTimeSlideUniformAgeConstrained.h"
#include "NodeTimeSlideUniformAgeConstrainedProposal.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlClade.h"
#include "RlTimeTree.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "ModelObject.h"
#include "Move.h"
#include "RbVector.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlMove.h"
#include "StochasticNode.h"

namespace RevBayesCore { class Clade; }
namespace RevBayesCore { class Proposal; }
namespace RevBayesCore { class Tree; }


using namespace RevLanguage;

Move_NodeTimeSlideUniformAgeConstrained::Move_NodeTimeSlideUniformAgeConstrained() : Move()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_NodeTimeSlideUniformAgeConstrained* Move_NodeTimeSlideUniformAgeConstrained::clone(void) const
{
    
    return new Move_NodeTimeSlideUniformAgeConstrained(*this);
}


void Move_NodeTimeSlideUniformAgeConstrained::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new move
    RevBayesCore::TypedDagNode<RevBayesCore::Tree> *tmp = static_cast<const TimeTree &>( tree->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::Tree> *t = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree> *>( tmp );
    
    const RevBayesCore::RbVector<RevBayesCore::Clade> &constr = static_cast<const ModelVector<Clade> &>( constraints->getRevObject() ).getValue();

    
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    size_t del = static_cast<const Natural &>( delay->getRevObject() ).getValue();
    
    RevBayesCore::Proposal *p = new RevBayesCore::NodeTimeSlideUniformAgeConstrainedProposal( t, constr );
    value = new RevBayesCore::MetropolisHastingsMove(p,w,del,false);
}


/** Get Rev type of object */
const std::string& Move_NodeTimeSlideUniformAgeConstrained::getClassType(void)
{
    
    static std::string rev_type = "Move_NodeTimeSlideUniformAgeConstrained";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& Move_NodeTimeSlideUniformAgeConstrained::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_NodeTimeSlideUniformAgeConstrained::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "NodeTimeSlideUniformAgeConstrained";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Move_NodeTimeSlideUniformAgeConstrained::getParameterRules(void) const
{
    
    static MemberRules member_rules;
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        
        member_rules.push_back( new ArgumentRule( "tree", TimeTree::getClassTypeSpec(), "The tree on which this move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        member_rules.push_back( new ArgumentRule( "constraints", ModelVector<Clade>::getClassTypeSpec(), "The age constraints.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );

        /* Inherit weight from Move, put it after variable */
        const MemberRules& inherited_rules = Move::getParameterRules();
        member_rules.insert( member_rules.end(), inherited_rules.begin(), inherited_rules.end() );
        
        rules_set = true;
    }
    
    return member_rules;
}

/** Get type spec */
const TypeSpec& Move_NodeTimeSlideUniformAgeConstrained::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}



/** Get type spec */
void Move_NodeTimeSlideUniformAgeConstrained::printValue(std::ostream &o) const
{
    
    o << "Move_NodeTimeSlideUniformAgeConstrained(";
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
void Move_NodeTimeSlideUniformAgeConstrained::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "tree" )
    {
        tree = var;
    }
    else if ( name == "constraints" )
    {
        constraints = var;
    }
    else
    {
        Move::setConstParameter(name, var);
    }
}
