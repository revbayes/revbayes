//
//  Func_featureInformedRates.cpp
//  revbayes-tensorphylo-proj
//
//  Created by Michael Landis on 7/12/22.
//  Copyright © 2022 Michael Landis. All rights reserved.
//

#include "Func_featureInformedRates.h"

#include <cstddef>

#include "FeatureInformedRateFunction.h"
#include "ModelVector.h"
#include "RlDeterministicNode.h"
#include "RlSimplex.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ModelObject.h"
#include "RbException.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "Simplex.h"
#include "StringUtilities.h"
#include "TypeSpec.h"

using namespace RevLanguage;

/** default constructor */
Func_featureInformedRates::Func_featureInformedRates( void ) : TypedFunction<ModelVector<ModelVector<RealPos> > >( ) {
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_featureInformedRates* Func_featureInformedRates::clone( void ) const {
    
    return new Func_featureInformedRates( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RbVector<RevBayesCore::RbVector<double> > >* Func_featureInformedRates::createFunction( void ) const
{
    
    // supplied arguments
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<RevBayesCore::RbVector<RevBayesCore::RbVector<long> > > >* cf = static_cast<const ModelVector<ModelVector<ModelVector<Natural> > >&>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<RevBayesCore::RbVector<RevBayesCore::RbVector<double> > > >* qf = static_cast<const ModelVector<ModelVector<ModelVector<Real> > >&>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* sigma = static_cast<const ModelVector<Real> &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* phi = static_cast<const ModelVector<Real> &>( this->args[3].getVariable()->getRevObject() ).getDagNode();
    
    // add error checking
    
    // create P matrix
    RevBayesCore::FeatureInformedRateFunction* f = new RevBayesCore::FeatureInformedRateFunction( cf, qf, sigma, phi );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_featureInformedRates::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "categoricalFeatures",
                                                  ModelVector<ModelVector<ModelVector<Natural> > >::getClassTypeSpec(),
                                                  "Vector of layers for categorical regional features.",
                                                  ArgumentRule::BY_CONSTANT_REFERENCE,
                                                  ArgumentRule::ANY ));
        argumentRules.push_back( new ArgumentRule( "quantitativeFeatures",
                                                  ModelVector<ModelVector<ModelVector<Real> > >::getClassTypeSpec(),
                                                  "Vector of layers for quantitative regional features.",
                                                  ArgumentRule::BY_CONSTANT_REFERENCE,
                                                  ArgumentRule::ANY ));
        argumentRules.push_back( new ArgumentRule( "sigma",
                                                  ModelVector<Real>::getClassTypeSpec(),
                                                  "Vector of effect parameters for each categorical feature layer.",
                                                  ArgumentRule::BY_CONSTANT_REFERENCE,
                                                  ArgumentRule::ANY ));
        argumentRules.push_back( new ArgumentRule( "phi",
                                                  ModelVector<Real>::getClassTypeSpec(),
                                                  "Vector of effect parameters for each quantitative feature layer.",
                                                  ArgumentRule::BY_CONSTANT_REFERENCE,
                                                  ArgumentRule::ANY ));
        
        rules_set = true;
    }
    
    return argumentRules;
}



const std::string& Func_featureInformedRates::getClassType(void)
{
    
    static std::string rev_type = "Func_featureInformedRates";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_featureInformedRates::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_featureInformedRates::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnFeatureInformedRates";
    
    return f_name;
}

const TypeSpec& Func_featureInformedRates::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
