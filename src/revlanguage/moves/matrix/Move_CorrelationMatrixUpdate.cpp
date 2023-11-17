/* 
 * File:   Move_CorrelationMatrixUpdate.cpp
 * Author: nl
 * 
 * Created on May 4, 2016
 */

#include "Move_CorrelationMatrixUpdate.h"

#include <cstddef>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "CorrelationMatrixProposal.h"
#include "RlBoolean.h"
#include "MetropolisHastingsMove.h"
#include "ModelVector.h"
#include "Real.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlMatrixReal.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "ModelObject.h"
#include "Move.h"
#include "RbBoolean.h"
#include "StochasticNode.h"
#include "StringUtilities.h"

namespace RevBayesCore { class MatrixReal; }
namespace RevBayesCore { class Proposal; }
namespace RevBayesCore { template <class valueType> class RbVector; }


using namespace RevLanguage;

Move_CorrelationMatrixUpdate::Move_CorrelationMatrixUpdate() : Move() {
    
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_CorrelationMatrixUpdate* Move_CorrelationMatrixUpdate::clone(void) const {
    
	return new Move_CorrelationMatrixUpdate(*this);
}


void Move_CorrelationMatrixUpdate::constructInternalObject( void ) {

    // we free the memory first
    delete value;
    
    // now allocate a new sliding move
    double l = static_cast<const RealPos &>( lambda->getRevObject() ).getValue();
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    RevBayesCore::TypedDagNode<RevBayesCore::MatrixReal >* tmp = static_cast<const MatrixReal &>( r->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::MatrixReal > *n = static_cast<RevBayesCore::StochasticNode<RevBayesCore::MatrixReal> *>( tmp );
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* sigma = static_cast<const ModelVector<RealPos> &>( s->getRevObject() ).getDagNode();
    bool t = static_cast<const RlBoolean &>( tune->getRevObject() ).getValue();
    
    RevBayesCore::Proposal *p = new RevBayesCore::CorrelationMatrixProposal(n, sigma, l);
    value = new RevBayesCore::MetropolisHastingsMove(p,w,t);

}


/** Get class name of object */
const std::string& Move_CorrelationMatrixUpdate::getClassType(void) {
    
    static std::string revClassType = "Move_CorrelationMatrixUpdate";
    
	return revClassType; 
}

/** Get class type spec describing type of object */
const TypeSpec& Move_CorrelationMatrixUpdate::getClassTypeSpec(void) {
    
    static TypeSpec revClassTypeSpec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
	return revClassTypeSpec; 
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_CorrelationMatrixUpdate::getMoveName( void ) const {

    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "CorrelationMatrixUpdate";
    
    return c_name;
}



/** Return member rules (no members) */
const MemberRules& Move_CorrelationMatrixUpdate::getParameterRules(void) const {
    
    static MemberRules move_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
        {
        move_member_rules.push_back( new ArgumentRule( "r"     , MatrixReal::getClassTypeSpec(),           "The correlation matrix on which this move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        move_member_rules.push_back( new ArgumentRule( "s"     , ModelVector<RealPos>::getClassTypeSpec(), "The vector of variances needed to complete the move.", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY ) );
        move_member_rules.push_back( new ArgumentRule( "lambda", RealPos::getClassTypeSpec(),              "The scaling factor (strength) of the proposal.", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new Real(1.0) ) );
        move_member_rules.push_back( new ArgumentRule( "tune"  , RlBoolean::getClassTypeSpec(),            "Should we tune the scaling factor during burnin?", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new RlBoolean( true ) ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
        }
    
    return move_member_rules;
}

/** Get type spec */
const TypeSpec& Move_CorrelationMatrixUpdate::getTypeSpec( void ) const {
    
    static TypeSpec type_spec = getClassTypeSpec();
    return type_spec;
}


/** Get type spec */
void Move_CorrelationMatrixUpdate::printValue(std::ostream &o) const {
    
    o << "Move_CorrelationMatrixUpdate(";
    if (r != NULL)
        {
        o << r->getName();
        }
    else
        {
        o << "?";
        }
    o << ")";
}


/** Set a member variable */
void Move_CorrelationMatrixUpdate::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
    if ( name == "r" )
        {
        r = var;
        }
    else if ( name == "s" )
        {
        s = var;
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
    else
        {
        Move::setConstParameter(name, var);
        }
}

