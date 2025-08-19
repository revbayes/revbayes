/* 
 * File:   Move_MatrixSingleElementScale.cpp
 * Author: nl
 * 
 * Created on 13 juillet 2014, 18:13
 */

#include "Move_MatrixSingleElementScale.h"

#include <cstddef>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RlBoolean.h"
#include "MatrixRealSingleElementScaleProposal.h"
#include "MetropolisHastingsMove.h"
#include "ModelVector.h"
#include "Real.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlMatrixReal.h"
#include "RlMatrixRealSymmetric.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "ModelObject.h"
#include "Move.h"
#include "RbBoolean.h"
#include "RbVector.h"
#include "StochasticNode.h"
#include "StringUtilities.h"

namespace RevBayesCore { class MatrixReal; }
namespace RevBayesCore { class Proposal; }


using namespace RevLanguage;

Move_MatrixSingleElementScale::Move_MatrixSingleElementScale() : Move()
{
    
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_MatrixSingleElementScale* Move_MatrixSingleElementScale::clone(void) const
{
    
	return new Move_MatrixSingleElementScale(*this);
}


void Move_MatrixSingleElementScale::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new sliding move
    double l = static_cast<const RealPos &>( lambda->getRevObject() ).getValue();
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    bool t = static_cast<const RlBoolean &>( tune->getRevObject() ).getValue();

    RevBayesCore::Proposal *p = NULL;

    if (v->getRevObject().isType( MatrixReal::getClassTypeSpec() ))
    {
        RevBayesCore::TypedDagNode<RevBayesCore::MatrixReal >* tmp = static_cast<const MatrixReal &>( v->getRevObject() ).getDagNode();
        RevBayesCore::StochasticNode<RevBayesCore::MatrixReal > *n = static_cast<RevBayesCore::StochasticNode<RevBayesCore::MatrixReal> *>( tmp );
        p = new RevBayesCore::MatrixRealSingleElementScaleProposal(n,l, v->getRevObject().isType( MatrixRealSymmetric::getClassTypeSpec() ) );
    }
    else if (v->getRevObject().isType( ModelVector<ModelVector<RealPos> >::getClassTypeSpec() ))
    {
        RevBayesCore::TypedDagNode<RevBayesCore::RbVector<RevBayesCore::RbVector<double> > >* tmp = static_cast<const ModelVector<ModelVector<RealPos> > &>( v->getRevObject() ).getDagNode();
        RevBayesCore::StochasticNode<RevBayesCore::RbVector<RevBayesCore::RbVector<double> > > *n = static_cast<RevBayesCore::StochasticNode<RevBayesCore::RbVector<RevBayesCore::RbVector<double> > > *>( tmp );
        p = new RevBayesCore::MatrixRealSingleElementScaleProposal(n,l);
    }
    else if (v->getRevObject().isType( ModelVector<ModelVector<Real> >::getClassTypeSpec() ))
    {
        RevBayesCore::TypedDagNode<RevBayesCore::RbVector<RevBayesCore::RbVector<double> > >* tmp = static_cast<const ModelVector<ModelVector<Real> > &>( v->getRevObject() ).getDagNode();
        RevBayesCore::StochasticNode<RevBayesCore::RbVector<RevBayesCore::RbVector<double> > > *n = static_cast<RevBayesCore::StochasticNode<RevBayesCore::RbVector<RevBayesCore::RbVector<double> > > *>( tmp );
        p = new RevBayesCore::MatrixRealSingleElementScaleProposal(n,l);
    }
    
    value = new RevBayesCore::MetropolisHastingsMove(p,w,t);

}


/** Get class name of object */
const std::string& Move_MatrixSingleElementScale::getClassType(void) {
    
    static std::string revClassType = "Move_MatrixElementScale";
    
	return revClassType; 
}

/** Get class type spec describing type of object */
const TypeSpec& Move_MatrixSingleElementScale::getClassTypeSpec(void)
{
    
    static TypeSpec revClassTypeSpec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
	return revClassTypeSpec; 
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_MatrixSingleElementScale::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "MatrixElementScale";
    
    return c_name;
}



/** Return member rules (no members) */
const MemberRules& Move_MatrixSingleElementScale::getParameterRules(void) const
{
    
    static MemberRules move_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        std::vector<TypeSpec> matTypes;
        matTypes.push_back( ModelVector<ModelVector<RealPos> >::getClassTypeSpec() );
        matTypes.push_back( ModelVector<ModelVector<Real> >::getClassTypeSpec() );
        matTypes.push_back( MatrixRealSymmetric::getClassTypeSpec() );
        matTypes.push_back( MatrixReal::getClassTypeSpec() );
        move_member_rules.push_back( new ArgumentRule( "x"     , matTypes, "The variable on which this move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        move_member_rules.push_back( new ArgumentRule( "lambda", RealPos::getClassTypeSpec()   , "The scaling factor (strength) of the proposal.", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new Real(1.0) ) );
        move_member_rules.push_back( new ArgumentRule( "tune"  , RlBoolean::getClassTypeSpec() , "Should we tune the scaling factor during burnin?", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new RlBoolean( true ) ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return move_member_rules;
}

/** Get type spec */
const TypeSpec& Move_MatrixSingleElementScale::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void Move_MatrixSingleElementScale::printValue(std::ostream &o) const
{
    
    o << "Move_MatrixElementScale(";
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
void Move_MatrixSingleElementScale::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
    if ( name == "x" )
    {
        v = var;
    }
    else if ( name == "lambda" )
    {
        lambda = var;
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

