#include "Func_PoMoTwo.h"
#include "Func_PoMoTwo.h"
#include "PoMoTwoRateMatrixFunction.h"
#include "ModelVector.h"
#include "Natural.h"
#include "RateMatrix_PoMoTwo.h"
#include "Real.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "RlRateMatrix.h"
#include "RlSimplex.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_PoMoTwo::Func_PoMoTwo( void ) : TypedFunction<RateMatrix>( ) {

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_PoMoTwo* Func_PoMoTwo::clone( void ) const {

    return new Func_PoMoTwo( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_PoMoTwo::createFunction( void ) const
{
    RevBayesCore::TypedDagNode<RevBayesCore::Simplex>* bf = static_cast<const Simplex &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* er = static_cast<const ModelVector<RealPos> &>( this->args[1].getVariable()->getRevObject() ).getDagNode();

    if ( er->getValue().size() != (bf->getValue().size() * (bf->getValue().size()-1) / 2.0) )
    {
        throw RbException("The dimensions between the base frequencies and the substitution rates do not match (they should be 4 and 6).");
    }

    RevBayesCore::TypedDagNode< long >* n = static_cast<const Natural &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    
    RevBayesCore::PoMoTwoRateMatrixFunction* f = new RevBayesCore::PoMoTwoRateMatrixFunction( n, er, bf );

    return f;
}


/* Get argument rules */
const ArgumentRules& Func_PoMoTwo::getArgumentRules( void ) const
{

  static ArgumentRules argumentRules = ArgumentRules();
  static bool          rules_set = false;

  if ( !rules_set )
  {
    argumentRules.push_back( new ArgumentRule( "baseFrequencies", Simplex::getClassTypeSpec(), "The stationary frequencies of the 4 nucleotide bases A, C, G and T.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    argumentRules.push_back( new ArgumentRule( "exchangeRates"      , ModelVector<Real>::getClassTypeSpec(), "The base exchangeabilities in this order AC, AG, AT, CG, CT and GT", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    argumentRules.push_back( new ArgumentRule( "Ne"    , Natural::getClassTypeSpec()          , "The effective population size", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

    rules_set = true;
  }

  return argumentRules;
}


const std::string& Func_PoMoTwo::getClassType(void)
{

  static std::string rev_type = "Func_PoMoTwo";

	return rev_type;

}


/* Get class type spec describing type of object */
const TypeSpec& Func_PoMoTwo::getClassTypeSpec(void)
{

  static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

	return rev_type_spec;

}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_PoMoTwo::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnPoMoTwo";

    return f_name;
}


const TypeSpec& Func_PoMoTwo::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
