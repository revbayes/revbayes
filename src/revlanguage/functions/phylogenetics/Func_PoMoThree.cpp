#include "Func_PoMoThree.h"
#include "Func_PoMoThree.h"
#include "PoMoThreeRateMatrixFunction.h"
#include "ModelVector.h"
#include "Natural.h"
#include "RateMatrix_PoMoThree.h"
#include "Real.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "RlRateMatrix.h"
#include "RlSimplex.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_PoMoThree::Func_PoMoThree( void ) : TypedFunction<RateMatrix>( ) {

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_PoMoThree* Func_PoMoThree::clone( void ) const {

    return new Func_PoMoThree( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_PoMoThree::createFunction( void ) const
{
    RevBayesCore::TypedDagNode<RevBayesCore::Simplex>* bf = static_cast<const Simplex &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* er = static_cast<const ModelVector<RealPos> &>( this->args[1].getVariable()->getRevObject() ).getDagNode();

    if ( er->getValue().size() != (bf->getValue().size() * (bf->getValue().size()-1) / 2.0) )
    {
        throw RbException("The dimensions between the base frequencies and the substitution rates do not match (they should be 4 and 6).");
    }

    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* sc = static_cast<const ModelVector<RealPos> &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< long >* n = static_cast<const Natural &>( this->args[3].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::PoMoThreeRateMatrixFunction* f = new RevBayesCore::PoMoThreeRateMatrixFunction( n, er, bf, sc );

    return f;
}


/* Get argument rules */
const ArgumentRules& Func_PoMoThree::getArgumentRules( void ) const
{

  static ArgumentRules argumentRules = ArgumentRules();
  static bool          rules_set = false;

  if ( !rules_set )
  {
    argumentRules.push_back( new ArgumentRule( "baseFrequencies", Simplex::getClassTypeSpec(), "The stationary frequencies of the 4 nucleotide bases A, C, G and T.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    argumentRules.push_back( new ArgumentRule( "exchangeRates"      , ModelVector<Real>::getClassTypeSpec(), "The base exchangeabilities in this order AC, AG, AT, CG, CT and GT", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    argumentRules.push_back( new ArgumentRule( "selectionCoefficients"      , ModelVector<Real>::getClassTypeSpec(), "The relative selection coefficients of C, G and T. A is set to 0.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    argumentRules.push_back( new ArgumentRule( "Ne"    , Natural::getClassTypeSpec()          , "The effective population size", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

    rules_set = true;
  }

  return argumentRules;
}


const std::string& Func_PoMoThree::getClassType(void)
{

    static std::string rev_type = "Func_PoMoThree";

	return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_PoMoThree::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
	return rev_type_spec;

}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_PoMoThree::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnPoMoThree";
    return f_name;
}


const TypeSpec& Func_PoMoThree::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();
    return type_spec;
}
