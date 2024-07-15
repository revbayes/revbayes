//
//  Func_PoMoReversibleMutationRates.cpp
//  RevBayesCore
//
//  Created by Rui Borges 
//  Copyright 2024
//

#include "Func_PoMoReversibleMutationRates.h"

#include "Natural.h"
#include "PoMoReversibleMutationRatesFunction.h"
#include "RlDeterministicNode.h"
#include "RlSimplex.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Real.h"
#include "RealPos.h"
#include "ModelVector.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "RlRateGenerator.h"
#include "TypeSpec.h"

namespace RevBayesCore { class RateGenerator; }

using namespace RevLanguage;

/** default constructor */
Func_PoMoReversibleMutationRates::Func_PoMoReversibleMutationRates( void ) : TypedFunction< ModelVector<RealPos> >( )
{
    
}


/** The clone function is a convenience function to create proper copies of inherited objected */
Func_PoMoReversibleMutationRates* Func_PoMoReversibleMutationRates::clone( void ) const
{
    
    return new Func_PoMoReversibleMutationRates( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RbVector<double> >* Func_PoMoReversibleMutationRates::createFunction( void ) const
{
    
    long na                                                             = static_cast<const Natural              &>( this->args[0].getVariable()->getRevObject() ).getValue();
    RevBayesCore::TypedDagNode<RevBayesCore::Simplex            >* pi   = static_cast<const Simplex              &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double>   >* rho  = static_cast<const ModelVector<RealPos> &>( this->args[2].getVariable()->getRevObject() ).getDagNode();

    if ( pi->getValue().size() != na )
    {
        throw RbException("The number of allele frequencies in pi does not match the number of alleles.");
    }
    if ( rho->getValue().size() != (na*na-na)/2 )
    {
        throw RbException("The number of alleles does not match the number of exchangeabilities: n_exchangeabilties = n_alleles*(n_alleles-1)/2");
    }

    RevBayesCore::PoMoReversibleMutationRatesFunction* pomomr = new RevBayesCore::PoMoReversibleMutationRatesFunction( na, pi, rho );
    
    return pomomr;
}


/* Get argument rules */
const ArgumentRules& Func_PoMoReversibleMutationRates::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set     = false;
    
    if ( !rules_set )
    {
        //Four arguments, root_base_frequencies, root_polymorphism_proportion, Q, virtual_population_size

        argumentRules.push_back( new ArgumentRule( "K"   , Natural::getClassTypeSpec(),               "Number of alleles", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "pi"  , Simplex::getClassTypeSpec(),               "Vector of base_frequencies: pi=[a0,a1,...,ak]", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "rho" , ModelVector<RealPos>::getClassTypeSpec(),  "Vector of exchangeabilities: rho=(rho_a0a1,rho_a0a2,...,akak-1)", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_PoMoReversibleMutationRates::getClassType(void)
{
    
    static std::string rev_type = "Func_PoMoReversibleMutationRates";
    
	return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_PoMoReversibleMutationRates::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_PoMoReversibleMutationRates::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnPoMoReversibleMutationRates";
    
    return f_name;
}


const TypeSpec& Func_PoMoReversibleMutationRates::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
