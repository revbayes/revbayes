//
//  Func_DECRates.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 3/3/15.
//  Copyright (c) 2015 Michael Landis. All rights reserved.
//

#include "DispersalExtinctionRateStructureFunction.h"
#include "Func_DECRates.h"
#include "Natural.h"
#include "RateMatrix_JC.h"
#include "Real.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "RlRateMatrix.h"
#include "RlSimplex.h"
#include "TransitionProbabilityMatrix.h"
#include "TypedDagNode.h"
#include <vector>

using namespace RevLanguage;

/** default constructor */
Func_DECRates::Func_DECRates( void ) : TypedFunction<ModelVector<RealPos> >( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_DECRates* Func_DECRates::clone( void ) const
{
    
    return new Func_DECRates( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RbVector<double> >* Func_DECRates::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RbVector<double> > >* dr = static_cast<const ModelVector<ModelVector<RealPos> > &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* er = static_cast<const ModelVector<RealPos> &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<long>* rs = static_cast<const Natural&>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    
    if (dr->getValue().size() != er->getValue().size())
    {
        RbException("ERROR: fnDispersalExtinction dispersalRates.size() != extinctionRates.size()!");
    }
    
    RevBayesCore::DispersalExtinctionRateStructureFunction* f = new RevBayesCore::DispersalExtinctionRateStructureFunction( dr,er,rs );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_DECRates::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "dispersalRates", ModelVector<ModelVector<RealPos> >::getClassTypeSpec(), "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "extinctionRates", ModelVector<RealPos>::getClassTypeSpec(), "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "maxRangeSize", Natural::getClassTypeSpec(), "", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(RbConstants::Integer::max) ) );

        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_DECRates::getClassType(void)
{
    
    static std::string rev_type = "Func_DECRates";
    
	return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_DECRates::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_DECRates::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnDECRates";
    
    return f_name;
}


const TypeSpec& Func_DECRates::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
