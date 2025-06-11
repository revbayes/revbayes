#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RlBoolean.h"
#include "MetropolisHastingsMove.h"
#include "ModelVector.h"
#include "Move_SynchronizedVectorFixedSingleElementSlide.h"
#include "Natural.h"
#include "RbException.h"
#include "Real.h"
#include "RealPos.h"
#include "RevObject.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "SynchronizedVectorFixedSingleElementSlideProposal.h"
#include "DagNode.h"
#include "ModelObject.h"
#include "Move.h"
#include "RbBoolean.h"
#include "RbVector.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlMove.h"
#include "StochasticNode.h"

namespace RevBayesCore { class Proposal; }


using namespace RevLanguage;

Move_SynchronizedVectorFixedSingleElementSlide::Move_SynchronizedVectorFixedSingleElementSlide() : Move()
{
    
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_SynchronizedVectorFixedSingleElementSlide* Move_SynchronizedVectorFixedSingleElementSlide::clone(void) const
{
    
    return new Move_SynchronizedVectorFixedSingleElementSlide(*this);
}


void Move_SynchronizedVectorFixedSingleElementSlide::constructInternalObject( void )
{
    
    // we free the memory first
    delete value;
    
    // now allocate a new sliding move
    double l = static_cast<const RealPos &>( lambda->getRevObject() ).getValue();
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    size_t del = static_cast<const Natural &>( delay->getRevObject() ).getValue();
    
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<RevBayesCore::RbVector<double> > >* tmp = static_cast<const ModelVector<ModelVector<RealPos> > &>( v->getRevObject() ).getDagNode();
    std::vector<const RevBayesCore::DagNode*> par = tmp->getParents();
    std::vector< RevBayesCore::StochasticNode<RevBayesCore::RbVector<double> > *> n;
    for (std::vector<const RevBayesCore::DagNode*>::const_iterator it = par.begin(); it != par.end(); ++it)
    {
        const RevBayesCore::StochasticNode<RevBayesCore::RbVector<double> > *the_node = dynamic_cast< const RevBayesCore::StochasticNode<RevBayesCore::RbVector<double> >* >( *it );
        if ( the_node != NULL )
        {
            n.push_back( const_cast< RevBayesCore::StochasticNode<RevBayesCore::RbVector<double> >* >( the_node ) );
        }
        else
        {
            throw RbException("Could not create a mvElementScale because the node isn't a vector of stochastic nodes.");
        }
    }
    
    
    bool t = static_cast<const RlBoolean &>( tune->getRevObject() ).getValue();
    
    std::vector<std::int64_t> e;
    if ( which_element->getRevObject().isType( ModelVector<Natural>::getClassTypeSpec() ) )
    {
        e = static_cast<const ModelVector<Natural> &>( which_element->getRevObject() ).getValue();
    }
    else
    {
        int index = (int)static_cast<const Natural &>( which_element->getRevObject() ).getValue();
        e.push_back( index );
    }
    
    // we need to offset the indices
    for (size_t i=0; i<e.size(); ++i)
    {
        --e[i];
    }

    
    RevBayesCore::Proposal *prop = new RevBayesCore::SynchronizedVectorFixedSingleElementSlideProposal(n, l, e);
    value = new RevBayesCore::MetropolisHastingsMove(prop, w, del, t);
    
}


/** Get Rev type of object */
const std::string& Move_SynchronizedVectorFixedSingleElementSlide::getClassType(void)
{
    
    static std::string rev_type = "Move_SynchronizedVectorFixedSingleElementSlide";
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& Move_SynchronizedVectorFixedSingleElementSlide::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_SynchronizedVectorFixedSingleElementSlide::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "SynchronizedVectorFixedSingleElementSlide";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Move_SynchronizedVectorFixedSingleElementSlide::getParameterRules(void) const
{
    
    static MemberRules move_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        move_member_rules.push_back( new ArgumentRule( "x"      , ModelVector<ModelVector<Real> >::getClassTypeSpec(), "The variable (a deterministic variable holding the vector of stochastic variable) on which this move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::DETERMINISTIC ) );
        
        std::vector<TypeSpec> index_types;
        index_types.push_back( Natural::getClassTypeSpec() );
        index_types.push_back( ModelVector<Natural>::getClassTypeSpec() );
        move_member_rules.push_back( new ArgumentRule( "element", index_types, "The index or indices of the element to scale.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        
        move_member_rules.push_back( new ArgumentRule( "lambda" , RealPos::getClassTypeSpec()          , "The scaling factor (strength) of this move.", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new Real(1.0) ) );
        move_member_rules.push_back( new ArgumentRule( "tune"   , RlBoolean::getClassTypeSpec()        , "Should we tune the scaling factor during burnin?", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new RlBoolean( true ) ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    return move_member_rules;
}

/** Get type spec */
const TypeSpec& Move_SynchronizedVectorFixedSingleElementSlide::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void Move_SynchronizedVectorFixedSingleElementSlide::printValue(std::ostream &o) const
{
    
    o << "Move_SynchronizedVectorFixedSingleElementSlide(";
    if (v != NULL)
    {
        o << v->getName();
    }
    else
    {
        o << "?";
    }
    o << ")";
    
}


/** 
 * Set a member variable
 */
void Move_SynchronizedVectorFixedSingleElementSlide::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    
    if ( name == "x" )
    {
        v = var;
    }
    else if ( name == "lambda" )
    {
        lambda = var;
    }
    else if ( name == "weight" )
    {
        weight = var;
    }
    else if ( name == "tune" )
    {
        tune = var;
    }
    else if ( name == "element" )
    {
        which_element = var;
    }
    else
    {
        Move::setConstParameter(name, var);
    }
    
}
