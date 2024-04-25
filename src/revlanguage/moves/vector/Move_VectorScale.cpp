#include <cstddef>
#include <ostream>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MetropolisHastingsMove.h"
#include "ModelVector.h"
#include "Move_VectorScale.h"
#include "RealPos.h"
#include "VectorScaleProposal.h"
#include "ModelObject.h"
#include "Move.h"
#include "RbBoolean.h"
#include "Real.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlBoolean.h"
#include "RlMove.h"
#include "StochasticNode.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class Proposal; }
namespace RevBayesCore { template <class valueType> class RbVector; }

using namespace RevLanguage;

Move_VectorScale::Move_VectorScale() : Move() {}


Move_VectorScale* Move_VectorScale::clone(void) const
{    
	return new Move_VectorScale(*this);
}


void Move_VectorScale::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new vector-scale move
    double l = static_cast<const RealPos &>( lambda->getRevObject() ).getValue();
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    size_t del = static_cast<const Natural &>( delay->getRevObject() ).getValue();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* tmp = static_cast<const ModelVector<RealPos> &>( x->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::RbVector<double> > *n = static_cast<RevBayesCore::StochasticNode<RevBayesCore::RbVector<double> > *>( tmp );
    
    bool t = static_cast<const RlBoolean &>( tune->getRevObject() ).getValue();
    
    RevBayesCore::Proposal *prop = new RevBayesCore::VectorScaleProposal(n,l);
    value = new RevBayesCore::MetropolisHastingsMove(prop,w,del,t);

}


const std::string& Move_VectorScale::getClassType(void)
{    
    static std::string rev_type = "Move_VectorScale";   
	return rev_type; 
}


const TypeSpec& Move_VectorScale::getClassTypeSpec(void)
{    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );    
	return rev_type_spec; 
}


std::string Move_VectorScale::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "VectorScale";    
    return c_name;
}


/**
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the VectorSingleElementScale move are:
 * (1) the variable which must be a positive real vector.
 * (2) the tuning parameter lambda that defines the size of the proposal (positive real)
 * (3) a flag whether auto-tuning should be used.
 *
 * \return The member rules.
 */
const MemberRules& Move_VectorScale::getParameterRules(void) const
{

    static MemberRules move_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        move_member_rules.push_back( new ArgumentRule( "x"     , ModelVector<RealPos>::getClassTypeSpec(), "The variable on which this move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        move_member_rules.push_back( new ArgumentRule( "lambda", RealPos::getClassTypeSpec()             , "The scaling parameter (strength) of the move.", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new Real(1.0) ) );
        move_member_rules.push_back( new ArgumentRule( "tune"  , RlBoolean::getClassTypeSpec()           , "Should we tune the scaling parameter during burnin?", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new RlBoolean( true ) ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return move_member_rules;
}


const TypeSpec& Move_VectorScale::getTypeSpec( void ) const {
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


void Move_VectorScale::printValue(std::ostream &o) const
{
    
    o << "Move_VectorScale(";
    if ( x != NULL )
    {
        o << x->getName();
    }
    else
    {
        o << "?";
    }
    o << ")";
}


void Move_VectorScale::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
    if ( name == "x" )
    {
        x = var;
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
