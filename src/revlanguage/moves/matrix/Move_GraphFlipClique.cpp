/*
 * File:   Move_GraphFlipClique.cpp
 * Author: nl
 *
 * Created on 13 juillet 2014, 18:13
 */

#include "Move_GraphFlipClique.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RlBoolean.h"
#include "GraphFlipCliqueProposal.h"
#include "MetropolisHastingsMove.h"
#include "ModelVector.h"
#include "Natural.h"
#include "Probability.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlMatrixRealSymmetric.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "Move.h"
#include "RbBoolean.h"
#include "RbVector.h"
#include "RlConstantNode.h"
#include "StochasticNode.h"
#include "StringUtilities.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

namespace RevBayesCore { class MatrixReal; }
namespace RevBayesCore { class Proposal; }


using namespace RevLanguage;

Move_GraphFlipClique::Move_GraphFlipClique() : Move()
{
    
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_GraphFlipClique* Move_GraphFlipClique::clone(void) const
{
    
    return new Move_GraphFlipClique(*this);
}


void Move_GraphFlipClique::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new sliding move
    double rs = static_cast<const RealPos&>( set_rate->getRevObject() ).getValue();
    double ep = static_cast<const Probability &>( edge_prob->getRevObject() ).getValue();
    double vp = static_cast<const Probability &>( vertex_prob->getRevObject() ).getValue();
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    bool t = static_cast<const RlBoolean &>( tune->getRevObject() ).getValue();
    const RevBayesCore::RbVector<std::int64_t>& e = static_cast<const ModelVector<Natural> &>( vertices->getRevObject() ).getValue();
    
    
    RevBayesCore::Proposal *p = NULL;
    
    RevBayesCore::TypedDagNode<RevBayesCore::MatrixReal >* tmp = static_cast<const MatrixRealSymmetric &>( v->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::MatrixReal > *n = static_cast<RevBayesCore::StochasticNode<RevBayesCore::MatrixReal> *>( tmp );
    
    bool symm = v->getRevObject().isType( MatrixRealSymmetric::getClassTypeSpec() );
    p = new RevBayesCore::GraphFlipCliqueProposal(n, e, rs, vp, ep, symm );
    
    value = new RevBayesCore::MetropolisHastingsMove(p,w,t);
    
}


/** Get class name of object */
const std::string& Move_GraphFlipClique::getClassType(void) {
    
    static std::string revClassType = "Move_GraphFlipClique";
    
    return revClassType;
}

/** Get class type spec describing type of object */
const TypeSpec& Move_GraphFlipClique::getClassTypeSpec(void)
{
    
    static TypeSpec revClassTypeSpec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
    return revClassTypeSpec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_GraphFlipClique::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "GraphFlipClique";
    
    return c_name;
}



/** Return member rules (no members) */
const MemberRules& Move_GraphFlipClique::getParameterRules(void) const
{
    
    static MemberRules move_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        std::vector<TypeSpec> matTypes;
        matTypes.push_back( MatrixRealSymmetric::getClassTypeSpec() );
        move_member_rules.push_back( new ArgumentRule( "x"     , matTypes, "The variable on which this move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        move_member_rules.push_back( new ArgumentRule( "r_set", RealPos::getClassTypeSpec(), "Rate that is equal to the expected number of vertex-sets to perturb (min 1).", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new RealPos(1.0) ) );
        move_member_rules.push_back( new ArgumentRule( "p_vertex", Probability::getClassTypeSpec(), "Probability of including a vertex in a vertex-set (min 2 vertices/set).", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new Probability(0.1) ) );
        move_member_rules.push_back( new ArgumentRule( "p_edge", Probability::getClassTypeSpec(), "Probability of forcing edge to consensus value (0,1) in a vertex-set.", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new Probability(1.0) ) );
        move_member_rules.push_back( new ArgumentRule( "vertices" , ModelVector<Natural>::getClassTypeSpec(), "A vector of vertices to target with this proposal. An empty vector is interpretted as the full list of vertices.", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY, new ModelVector<Natural>() ) );
        move_member_rules.push_back( new ArgumentRule( "tune"  , RlBoolean::getClassTypeSpec() , "Should we tune the scaling factor during burnin?", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new RlBoolean( true ) ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return move_member_rules;
}

/** Get type spec */
const TypeSpec& Move_GraphFlipClique::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void Move_GraphFlipClique::printValue(std::ostream &o) const
{
    
    o << "Move_GraphFlipClique(";
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
void Move_GraphFlipClique::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
    if ( name == "x" )
    {
        v = var;
    }
    else if ( name == "vertices" )
    {
        vertices = var;
    }
    else if ( name == "r_set" )
    {
        set_rate = var;
    }
    else if ( name == "p_vertex" )
    {
        vertex_prob = var;
    }
    else if ( name == "p_edge" )
    {
        edge_prob = var;
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

