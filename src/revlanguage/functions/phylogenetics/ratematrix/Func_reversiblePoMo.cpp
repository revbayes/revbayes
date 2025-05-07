#include "Func_reversiblePoMo.h"
//
//  Func_reversiblePoMo.cpp
//  RevBayesCore
//
//  Created by Bastien Boussau on 22/8/14.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include "ReversiblePoMoRateMatrixFunction.h"
#include "ModelVector.h"
#include "Natural.h"
#include "Real.h"
#include "RlDeterministicNode.h"
#include "RlRateMatrix.h"
#include "RlSimplex.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ModelObject.h"
#include "RbException.h"
#include "RbVector.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "Simplex.h"
#include "TypeSpec.h"

namespace RevLanguage { class RealPos; }

using namespace RevLanguage;

/** default constructor */
Func_reversiblePoMo::Func_reversiblePoMo( void ) : TypedFunction<RateMatrix>( ) {

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_reversiblePoMo* Func_reversiblePoMo::clone( void ) const {

    return new Func_reversiblePoMo( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_reversiblePoMo::createFunction( void ) const
{
    RevBayesCore::TypedDagNode<RevBayesCore::Simplex>* bf = static_cast<const Simplex &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* er = static_cast<const ModelVector<RealPos> &>( this->args[1].getVariable()->getRevObject() ).getDagNode();

    if ( er->getValue().size() != (bf->getValue().size() * (bf->getValue().size()-1) / 2.0) )
    {
        throw RbException("The dimensions between the base frequencies and the substitution rates do not match (they should be 4 and 6).");
    }


  //RevBayesCore::TypedDagNode<RevBayesCore::RateGenerator >* q = static_cast<const RateMatrix &>( this->args[0].getVariable()->getRevObject() ).getDagNode();

    RevBayesCore::TypedDagNode< std::int64_t >* n = static_cast<const Natural &>( this->args[2].getVariable()->getRevObject() ).getDagNode();

    RevBayesCore::ReversiblePoMoRateMatrixFunction* f = new RevBayesCore::ReversiblePoMoRateMatrixFunction( n, er, bf );

    return f;
}


/* Get argument rules */
const ArgumentRules& Func_reversiblePoMo::getArgumentRules( void ) const
{

  static ArgumentRules argumentRules = ArgumentRules();
  static bool          rules_set = false;

  if ( !rules_set )
  {
    argumentRules.push_back( new ArgumentRule( "baseFrequencies", Simplex::getClassTypeSpec(), "The stationary frequencies of the 4 DNA states.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    argumentRules.push_back( new ArgumentRule( "exchangeRates"      , ModelVector<Real>::getClassTypeSpec(), "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

    //argumentRules.push_back( new ArgumentRule( "mutationRates", RateGenerator::getClassTypeSpec()    , "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    argumentRules.push_back( new ArgumentRule( "virtualNe"    , Natural::getClassTypeSpec()          , "The virtual population size", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

    rules_set = true;
  }

  return argumentRules;
}


const std::string& Func_reversiblePoMo::getClassType(void)
{

    static std::string rev_type = "Func_reversiblePoMo";

	return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_reversiblePoMo::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_reversiblePoMo::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnReversiblePoMo";

    return f_name;
}


const TypeSpec& Func_reversiblePoMo::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
