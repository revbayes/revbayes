#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MetropolisHastingsMove.h"
#include "ModelVector.h"
#include "Move_VectorSlide.h"
#include "Natural.h"
#include "Real.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlBoolean.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "VectorSlideProposal.h"
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

Move_VectorSlide::Move_VectorSlide() : Move()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_VectorSlide* Move_VectorSlide::clone(void) const
{
    
	return new Move_VectorSlide(*this);
}


void Move_VectorSlide::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new vector-slide move
    const RevBayesCore::RbVector<std::int64_t> &e = static_cast<const ModelVector<Natural> &>( inidices->getRevObject() ).getValue();
    double l = static_cast<const RealPos &>( delta->getRevObject() ).getValue();
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* tmp = static_cast<const ModelVector<RealPos> &>( x->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::RbVector<double> > *n = static_cast<RevBayesCore::StochasticNode<RevBayesCore::RbVector<double> > *>( tmp );
    
    bool t = static_cast<const RlBoolean &>( tune->getRevObject() ).getValue();
    
    RevBayesCore::Proposal *prop = new RevBayesCore::VectorSlideProposal(n,e,l);
    value = new RevBayesCore::MetropolisHastingsMove(prop,w,t);
}


/** Get Rev type of object */
const std::string& Move_VectorSlide::getClassType(void)
{
    
    static std::string rev_type = "Move_VectorSlide";
    
	return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& Move_VectorSlide::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_VectorSlide::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "VectorSlide";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Move_VectorSlide::getParameterRules(void) const
{
    
    static MemberRules move_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        move_member_rules.push_back( new ArgumentRule( "x"        , ModelVector<Real>::getClassTypeSpec()   , "The variable on which this move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        move_member_rules.push_back( new ArgumentRule( "elements" , ModelVector<Natural>::getClassTypeSpec(), "The indices of elements. All are used if this is empty.", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY, new ModelVector<Natural>() ) );
        move_member_rules.push_back( new ArgumentRule( "delta"    , RealPos::getClassTypeSpec()             , "The window size parameter.", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new RealPos(1.0) ) );
        move_member_rules.push_back( new ArgumentRule( "tune"     , RlBoolean::getClassTypeSpec()           , "Should we tune the window size during burnin?", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new RlBoolean( true ) ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return move_member_rules;
}

/** Get type spec */
const TypeSpec& Move_VectorSlide::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void Move_VectorSlide::printValue(std::ostream &o) const
{
    
    o << "Move_VectorSlide(";
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
void Move_VectorSlide::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
    if ( name == "x" )
    {
        x = var;
    }
    else if ( name == "elements" )
    {
        inidices = var;
    }
    else if ( name == "delta" )
    {
        delta = var;
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
