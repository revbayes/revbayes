#include "Func_revPoMoThree4.h"
#include "Func_revPoMoThree4.h"
#include "revPoMoThree4RateMatrixFunction.h"
#include "ModelVector.h"
#include "Natural.h"
#include "RateMatrix_revPoMoThree4.h"
#include "Real.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "RlRateMatrix.h"
#include "RlSimplex.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_revPoMoThree4::Func_revPoMoThree4( void ) : TypedFunction<RateMatrix>( ) {

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_revPoMoThree4* Func_revPoMoThree4::clone( void ) const {

    return new Func_revPoMoThree4( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_revPoMoThree4::createFunction( void ) const
{
    RevBayesCore::TypedDagNode<RevBayesCore::Simplex>*           bf = static_cast<const Simplex              &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* ex = static_cast<const ModelVector<RealPos> &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* fc = static_cast<const ModelVector<RealPos> &>( this->args[2].getVariable()->getRevObject() ).getDagNode();

    if ( bf->getValue().size() !=  4 )
    {
        throw RbException("The number of base frequencies pi does not match the number of alleles: 4.");
    }

    if ( ex->getValue().size() !=  6 )
    {
        throw RbException("The number of exchangeabilities rho does not match the number of pairwise combinations of alleles: 6.");
    }
    
    if ( fc->getValue().size() !=  4 )
    {
        throw RbException("The number of fitness coefficients phi does not match the number of alleles: 4.");
    }

    RevBayesCore::revPoMoThree4RateMatrixFunction* f = new RevBayesCore::revPoMoThree4RateMatrixFunction( bf, ex, fc );

    return f;
}


/* Get argument rules */
const ArgumentRules& Func_revPoMoThree4::getArgumentRules( void ) const
{

  static ArgumentRules argumentRules = ArgumentRules();
  static bool          rules_set = false;

  if ( !rules_set )
  {

    argumentRules.push_back( new ArgumentRule( "pi"         , Simplex::getClassTypeSpec(), "Vector of allele frequencies: pi=(pi_a0,pi_a1,pi_a2,pi_a3)", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    argumentRules.push_back( new ArgumentRule( "rho"        , ModelVector<RealPos>::getClassTypeSpec(), "Vector of exchangeabilities: rho=(rho_a0a1,rho_a0a2,rho_a0a3,rho_a1a2,rho_a1a3,rho_a2a3)", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    argumentRules.push_back( new ArgumentRule( "phi"        , ModelVector<RealPos>::getClassTypeSpec(), "Vector of fitness coefficients: phi=(phi_a0,phi_a1,phi_a2,phi_a3)", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
 
    rules_set = true;
  }

  return argumentRules;
}


const std::string& Func_revPoMoThree4::getClassType(void)
{

  static std::string rev_type = "Func_revPoMoThree4";

	return rev_type;

}


/* Get class type spec describing type of object */
const TypeSpec& Func_revPoMoThree4::getClassTypeSpec(void)
{

  static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

	return rev_type_spec;

}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_revPoMoThree4::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnReversiblePoMoThree4";

    return f_name;
}


const TypeSpec& Func_revPoMoThree4::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
