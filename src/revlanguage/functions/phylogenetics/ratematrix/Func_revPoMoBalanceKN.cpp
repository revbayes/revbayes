#include "Func_revPoMoBalanceKN.h"
#include "revPoMoBalanceKNRateMatrixFunction.h"
#include "RateMatrix_revPoMoBalanceKN.h"
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
Func_revPoMoBalanceKN::Func_revPoMoBalanceKN( void ) : TypedFunction<RateMatrix>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_revPoMoBalanceKN* Func_revPoMoBalanceKN::clone( void ) const
{
    
    return new Func_revPoMoBalanceKN( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_revPoMoBalanceKN::createFunction( void ) const
{
    RevBayesCore::TypedDagNode< long                          >* na = static_cast<const Natural              &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< long                          >* ni = static_cast<const Natural              &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::Simplex          >* p  = static_cast<const Simplex              &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* r  = static_cast<const ModelVector<RealPos> &>( this->args[3].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* s  = static_cast<const ModelVector<Real>    &>( this->args[4].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* b  = static_cast<const ModelVector<RealPos> &>( this->args[5].getVariable()->getRevObject() ).getDagNode();

    if ( ni->getValue() % 2 != 0 )
    {
        throw RbException("The population size in the reversible balancing selection model must be even.");
    }

    if ( p->getValue().size() !=  na->getValue() )
    {
        throw RbException("The number of base frequencies pi does not match the number of alleles: n_base_freq = n_alleles.");
    }

    if ( r->getValue().size() != ( (na->getValue())*(na->getValue())-na->getValue() )/2 )
    {
        throw RbException("The number of exchangeabilities rho does not match the number of pairwise combinations of alleles: n_exchangeabilities = (n_alleles*n_alleles-n_alleles)/2.");
    }

    if ( s->getValue().size() !=  na->getValue() )
    {
        throw RbException("The number of fitness coefficients phi does not match the number of alleles: n_fit_coeff = n_alleles.");
    }

    if ( b->getValue().size() != ( (na->getValue())*(na->getValue())-na->getValue() )/2 )
    {
        throw RbException("The number of balancing selection coefficients beta does not match the number of pairwise combinations of alleles: n_bal_sel_str = (n_alleles*n_alleles-n_alleles)/2.");
    }

    RevBayesCore::revPoMoBalanceKNRateMatrixFunction* f = new RevBayesCore::revPoMoBalanceKNRateMatrixFunction( na, ni, p, r, s, b );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_revPoMoBalanceKN::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "K"          , Natural::getClassTypeSpec(), "Number of alleles", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "N"          , Natural::getClassTypeSpec(), "Population size", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "pi"         , Simplex::getClassTypeSpec(), "Vector of allele frequencies: pi=(pi_a0,pi_a1,...,pi_ak)", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "rho"        , ModelVector<RealPos>::getClassTypeSpec(), "Vector of exchangeabilities: rho=(rho_a0a1,rho_a0a2,...,rho_a0ak,rho_a1a2,...)", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "phi"        , ModelVector<Real>   ::getClassTypeSpec(), "Vector of fitness coefficients: phi=(phi_0,phi_1,...,phi_ak)", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "beta"       , ModelVector<RealPos>::getClassTypeSpec(), "Vector of balancing selection coefficients: beta=(beta_a0a1,beta_a0a2,...,beta_a0ak,beta_a1a2,...)", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_revPoMoBalanceKN::getClassType(void)
{
    
    static std::string rev_type = "Func_revPoMoBalanceKN";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_revPoMoBalanceKN::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_revPoMoBalanceKN::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnReversiblePoMoBalanceKN";
    
    return f_name;
}


const TypeSpec& Func_revPoMoBalanceKN::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}



