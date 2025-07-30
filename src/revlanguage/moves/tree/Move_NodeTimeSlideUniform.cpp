#include <cstddef>
#include <ostream>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MetropolisHastingsMove.h"
#include "ModelVector.h"
#include "Move_NodeTimeSlideUniform.h"
#include "NodeTimeSlideUniformProposal.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlTimeTree.h"
#include "TypeSpec.h"
#include "Move.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlMove.h"
#include "StochasticNode.h"


using namespace RevLanguage;

Move_NodeTimeSlideUniform::Move_NodeTimeSlideUniform() : Move()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_NodeTimeSlideUniform* Move_NodeTimeSlideUniform::clone(void) const
{
    
	return new Move_NodeTimeSlideUniform(*this);
}


void Move_NodeTimeSlideUniform::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new move
    
    // get the tree(s) variable
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
        throw RbException() << "Wrong tree type '" << tree->getRevObject().getType() << "'.";
    }
    
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    size_t del = static_cast<const Natural &>( delay->getRevObject() ).getValue();
    
    RevBayesCore::Proposal *p = new RevBayesCore::NodeTimeSlideUniformProposal( t, vec_t );
    value = new RevBayesCore::MetropolisHastingsMove(p,w,del,false);
}


/** Get Rev type of object */
const std::string& Move_NodeTimeSlideUniform::getClassType(void)
{
    
    static std::string rev_type = "Move_NodeTimeSlideUniform";
    
	return rev_type; 
}

/** Get class type spec describing type of object */
const TypeSpec& Move_NodeTimeSlideUniform::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_NodeTimeSlideUniform::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "NodeTimeSlideUniform";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Move_NodeTimeSlideUniform::getParameterRules(void) const
{
    
    static MemberRules member_rules;
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        std::vector<TypeSpec> tree_var_types;
        tree_var_types.push_back( TimeTree::getClassTypeSpec() );
        tree_var_types.push_back( ModelVector<TimeTree>::getClassTypeSpec() );
        
        member_rules.push_back( new ArgumentRule( "tree", tree_var_types, "The tree on which this move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inherited_rules = Move::getParameterRules();
        member_rules.insert( member_rules.end(), inherited_rules.begin(), inherited_rules.end() );
        
        rules_set = true;
    }
    
    return member_rules;
}

/** Get type spec */
const TypeSpec& Move_NodeTimeSlideUniform::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}



/** Get type spec */
void Move_NodeTimeSlideUniform::printValue(std::ostream &o) const
{
    
    o << "Move_NodeTimeSlideUniform(";
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
void Move_NodeTimeSlideUniform::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "tree" )
    {
        tree = var;
    }
    else
    {
        Move::setConstParameter(name, var);
    }
}
