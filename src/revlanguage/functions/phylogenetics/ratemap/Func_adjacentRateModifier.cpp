//
//  Func_adjacentRateModifier.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 2/3/17.
//  Copyright © 2017 Michael Landis. All rights reserved.
//

#include "Func_adjacentRateModifier.h"

#include <cstddef>
#include <cstdint>
#include <ostream>

#include "AdjacentRateModifierFunction.h"
#include "CharacterHistoryRateModifier.h"
#include "MatrixReal.h"
#include "ModelVector.h"
#include "RbVector.h"
#include "Real.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "RlMatrixReal.h"
#include "RlMatrixRealSymmetric.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ModelObject.h"
#include "Natural.h"
#include "RbException.h"
#include "RevObject.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "StochasticNode.h"
#include "TypeSpec.h"

using namespace RevLanguage;

/** default constructor */
Func_adjacentRateModifier::Func_adjacentRateModifier( void ) : TypedFunction< CharacterHistoryRateModifier >( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_adjacentRateModifier* Func_adjacentRateModifier::clone( void ) const {
    
    return new Func_adjacentRateModifier( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::CharacterHistoryRateModifier >* Func_adjacentRateModifier::createFunction( void ) const
{
    
    RevPtr<const RevVariable> gainFactor = this->args[0].getVariable();
    RevPtr<const RevVariable> lossFactor = this->args[1].getVariable();
    RevPtr<const RevVariable> width      = this->args[2].getVariable();
    RevPtr<const RevVariable> matrix     = this->args[3].getVariable();
    RevPtr<const RevVariable> numStates  = this->args[4].getVariable();
    RevPtr<const RevVariable> numChars   = this->args[5].getVariable();
    
    RevBayesCore::TypedDagNode<double>* gf = static_cast<const Real&>( gainFactor->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* lf = static_cast<const Real&>( lossFactor->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<std::int64_t>* w = static_cast<const Natural&>( width->getRevObject() ).getDagNode();
    
    size_t ns = static_cast<const Natural&>( numStates->getRevObject() ).getValue();
    size_t nc = static_cast<const Natural&>( numChars->getRevObject() ).getValue();

//    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<RevBayesCore::RbVector<std::int64_t> > >* c = NULL;
//    if ( this->args[3].getVariable()->getRevObject().isType( ModelVector<ModelVector<Natural> >::getClassTypeSpec() ) )
//    {
//        
//        c = static_cast<const ModelVector<ModelVector<Natural> > &>( this->args[3].getVariable()->getRevObject() ).getDagNode();
//    }
    
    RevBayesCore::AdjacentRateModifierFunction* f = NULL;
    
    if (matrix->getRevObject().isType( MatrixReal::getClassTypeSpec() ))
    {
        RevBayesCore::TypedDagNode<RevBayesCore::MatrixReal >* tmp = static_cast<const MatrixReal &>( matrix->getRevObject() ).getDagNode();
        RevBayesCore::StochasticNode<RevBayesCore::MatrixReal > *n = static_cast<RevBayesCore::StochasticNode<RevBayesCore::MatrixReal> *>( tmp );
        
        // validate input
        if (n != NULL) {
            size_t d1 = n->getValue().size();
            size_t d2 = n->getValue()[0].size();
            if (d1 != nc || d2 != nc)
            {
                std::stringstream ss;
                ss << "Context matrix dimensions (" << d1 << " x " << d2 << ") do not match number of characters (" << nc << ")";
                throw RbException(ss.str());
            }
        }

        f = new RevBayesCore::AdjacentRateModifierFunction(gf, lf, w, n, ns, nc);
//        p = new RevBayesCore::MatrixRealSingleElementSlideProposal(n,l, v->getRevObject().isType( MatrixRealSymmetric::getClassTypeSpec() ) );
    }
    else
    {
        throw RbException("Support for Real[][] coming soon.");
    }
//    else if (matrix->getRevObject().isType( ModelVector<ModelVector<RealPos> >::getClassTypeSpec() ))
//    {
//        RevBayesCore::TypedDagNode<RevBayesCore::RbVector<RevBayesCore::RbVector<double> > >* tmp = static_cast<const ModelVector<ModelVector<RealPos> > &>( v->getRevObject() ).getDagNode();
//        RevBayesCore::StochasticNode<RevBayesCore::RbVector<RevBayesCore::RbVector<double> > > *n = static_cast<RevBayesCore::StochasticNode<RevBayesCore::RbVector<RevBayesCore::RbVector<double> > > *>( tmp );
//        p = new RevBayesCore::MatrixRealSingleElementSlideProposal(n,l);
//    }
//    else if (matrix->getRevObject().isType( ModelVector<ModelVector<Real> >::getClassTypeSpec() ))
//    {
//        RevBayesCore::TypedDagNode<RevBayesCore::RbVector<RevBayesCore::RbVector<double> > >* tmp = static_cast<const ModelVector<ModelVector<Real> > &>( v->getRevObject() ).getDagNode();
//        RevBayesCore::StochasticNode<RevBayesCore::RbVector<RevBayesCore::RbVector<double> > > *n = static_cast<RevBayesCore::StochasticNode<RevBayesCore::RbVector<RevBayesCore::RbVector<double> > > *>( tmp );
//        p = new RevBayesCore::MatrixRealSingleElementSlideProposal(n,l);
//    }
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_adjacentRateModifier::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rulesSet = false;
    
    if ( !rulesSet )
    {
        
        std::vector<TypeSpec> matTypes;
        matTypes.push_back( ModelVector<ModelVector<RealPos> >::getClassTypeSpec() );
        matTypes.push_back( ModelVector<ModelVector<Real> >::getClassTypeSpec() );
        matTypes.push_back( MatrixRealSymmetric::getClassTypeSpec() );
        matTypes.push_back( MatrixReal::getClassTypeSpec() );
//        move_member_rules.push_back( new ArgumentRule( "x"     , matTypes, "The variable on which this move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        argumentRules.push_back( new ArgumentRule( "gainFactor", Real::getClassTypeSpec(), "Multiplicative factor (r' = r * e^{ n_1 * f }) for characters in context set", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Real(0.0) ) );
        argumentRules.push_back( new ArgumentRule( "lossFactor", Real::getClassTypeSpec(), "Multiplicative factor (r' = r * e^{ n_0 * f }) for characters in context set", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Real(0.0) ) );
        argumentRules.push_back( new ArgumentRule( "width",  Natural::getClassTypeSpec(), "Width of context-dependence window", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Natural(1) ) );
        argumentRules.push_back( new ArgumentRule( "matrix",  matTypes, "Weighted character adjacency matrix.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        argumentRules.push_back( new ArgumentRule( "numStates", Natural::getClassTypeSpec(), "Number of states", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY) );
        argumentRules.push_back( new ArgumentRule( "numChars", Natural::getClassTypeSpec(), "Number of characters", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY) );
        
        rulesSet = true;
    }
    
    return argumentRules;
}


const std::string& Func_adjacentRateModifier::getClassType(void)
{
    
    static std::string revType = "Func_adjacentRateModifier";
    
    return revType;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_adjacentRateModifier::getClassTypeSpec(void)
{
    
    static TypeSpec revTypeSpec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return revTypeSpec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_adjacentRateModifier::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnAdjacentRateModifier";
    
    return f_name;
}


const TypeSpec& Func_adjacentRateModifier::getTypeSpec( void ) const
{
    
    static TypeSpec typeSpec = getClassTypeSpec();
    
    return typeSpec;
}

