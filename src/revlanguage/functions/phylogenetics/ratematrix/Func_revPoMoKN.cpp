#include "Func_revPoMoKN.h"
#include "revPoMoKNRateMatrixFunction.h"
#include "RateMatrix_revPoMoKN.h"
#include "ModelVector.h"
#include "Natural.h"
#include "Real.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "RlRateMatrix.h"
#include "RlSimplex.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_revPoMoKN::Func_revPoMoKN( void ) : TypedFunction<RateMatrix>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_revPoMoKN* Func_revPoMoKN::clone( void ) const
{
    
    return new Func_revPoMoKN( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_revPoMoKN::createFunction( void ) const
{
    RevBayesCore::TypedDagNode< std::int64_t                          >* na = static_cast<const Natural              &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< std::int64_t                          >* ni = static_cast<const Natural              &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::Simplex          >* bf = static_cast<const Simplex              &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* ex = static_cast<const ModelVector<RealPos> &>( this->args[3].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* fc = static_cast<const ModelVector<RealPos> &>( this->args[4].getVariable()->getRevObject() ).getDagNode();

    if ( ex->getValue().size() != ( (na->getValue())*(na->getValue())-na->getValue() )/2 )
    {
        throw RbException("The number of alleles does not match the number of exchangeabilities given: n_exchangeabilities = (n_alleles*n_alleles-n_alleles)/2");
    }
    if ( bf->getValue().size() != na->getValue() )
    {
        throw RbException("The number of alleles does not match the number of allele frequencies given: n_allele_freq = n_alleles");
    }
    if ( fc->getValue().size() != na->getValue() )
    {
        throw RbException("The number of alleles does not match the number of fitness coefficients given: n_fit_coeff = n_alleles");
    }

    RevBayesCore::revPoMoKNRateMatrixFunction* f = new RevBayesCore::revPoMoKNRateMatrixFunction( na, ni, bf, ex, fc );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_revPoMoKN::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "K"   , Natural::getClassTypeSpec(), "Number of alleles", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "N"   , Natural::getClassTypeSpec(), "Number of individuals", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "pi"  , Simplex::getClassTypeSpec(), "Vector of allele frequencies_ pi=(pi_a0,pi_a1,...,pi_ak)", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "rho" , ModelVector<RealPos>::getClassTypeSpec(), "Vector of exchangeabilities: rho=(rho_a0a1,rho_a0a2,...,rho_a0ak,rho_a1a2,...)", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "phi" , ModelVector<RealPos>::getClassTypeSpec(), "Vector of fitness coefficients: phi=(phi_0,phi_1,...,phi_ak)", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );  

        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_revPoMoKN::getClassType(void)
{
    
    static std::string rev_type = "Func_revPoMoKN";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_revPoMoKN::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_revPoMoKN::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnReversiblePoMoKN";
    
    return f_name;
}


const TypeSpec& Func_revPoMoKN::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


