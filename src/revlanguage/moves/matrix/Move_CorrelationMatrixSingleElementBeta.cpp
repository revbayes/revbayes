/* 
 * File:   Move_CorrelationMatrixSingleElementBeta.cpp
 * Author: Michael R. May
 *
 * Created on 5 August 2017
 */

#include "Move_CorrelationMatrixSingleElementBeta.h"

#include <cstddef>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RlBoolean.h"
#include "CorrelationMatrixElementBetaProposal.h"
#include "MetropolisHastingsMove.h"
#include "Probability.h"
#include "Real.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlMatrixReal.h"
#include "TypeSpec.h"
#include "Move.h"
#include "RbBoolean.h"
#include "RbHelpReference.h"
#include "StochasticNode.h"
#include "StringUtilities.h"

namespace RevBayesCore { class MatrixReal; }
namespace RevBayesCore { class Proposal; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }


using namespace RevLanguage;

Move_CorrelationMatrixSingleElementBeta::Move_CorrelationMatrixSingleElementBeta() : Move()
{
    
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_CorrelationMatrixSingleElementBeta* Move_CorrelationMatrixSingleElementBeta::clone(void) const
{
    
	return new Move_CorrelationMatrixSingleElementBeta(*this);
}


void Move_CorrelationMatrixSingleElementBeta::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new sliding move
    double a = static_cast<const RealPos &>( alpha->getRevObject() ).getValue();
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    double r = static_cast<const Probability &>( tuneTarget->getRevObject() ).getValue();
    RevBayesCore::TypedDagNode<RevBayesCore::MatrixReal >* tmp = static_cast<const MatrixReal &>( v->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::MatrixReal > *n = static_cast<RevBayesCore::StochasticNode<RevBayesCore::MatrixReal> *>( tmp );
    bool t = static_cast<const RlBoolean &>( tune->getRevObject() ).getValue();
    
    RevBayesCore::Proposal *p = new RevBayesCore::CorrelationMatrixElementBetaProposal(n,a,r);
    value = new RevBayesCore::MetropolisHastingsMove(p,w,t);

}


/** Get class name of object */
const std::string& Move_CorrelationMatrixSingleElementBeta::getClassType(void) { 
    
    static std::string revClassType = "Move_CorrelationMatrixSingleElementBeta";
    
	return revClassType; 
}

/** Get class type spec describing type of object */
const TypeSpec& Move_CorrelationMatrixSingleElementBeta::getClassTypeSpec(void)
{
    
    static TypeSpec revClassTypeSpec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
	return revClassTypeSpec; 
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_CorrelationMatrixSingleElementBeta::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "CorrelationMatrixSingleElementBeta";
    
    return c_name;
}



/** Return member rules (no members) */
const MemberRules& Move_CorrelationMatrixSingleElementBeta::getParameterRules(void) const
{
    
    static MemberRules move_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        move_member_rules.push_back( new ArgumentRule( "x"     , MatrixReal::getClassTypeSpec(), "The variable on which this move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        move_member_rules.push_back( new ArgumentRule( "alpha",  RealPos::getClassTypeSpec()   , "The shape of the beta distribution to draw from.", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new Real(10.0) ) );
        move_member_rules.push_back( new ArgumentRule( "tune"  , RlBoolean::getClassTypeSpec() , "Should we tune the scaling factor during burnin?", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new RlBoolean( true ) ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return move_member_rules;
}

/** Get type spec */
const TypeSpec& Move_CorrelationMatrixSingleElementBeta::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void Move_CorrelationMatrixSingleElementBeta::printValue(std::ostream &o) const
{
    
    o << "Move_CorrelationMatrixSingleElementBeta(";
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


/** Set a member variable */
void Move_CorrelationMatrixSingleElementBeta::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
    if ( name == "x" )
    {
        v = var;
    }
    else if ( name == "alpha" )
    {
        alpha = var;
    }
    else if ( name == "weight" ) {
        weight = var;
    }
    else if ( name == "tune" ) {
        tune = var;
    }
    else {
        Move::setConstParameter(name, var);
    }
}

