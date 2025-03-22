#include "Func_PoMoKN.h"
#include "PoMoKNRateMatrixFunction.h"
#include "RateMatrix_PoMoKN.h"
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
Func_PoMoKN::Func_PoMoKN( void ) : TypedFunction<RateMatrix>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_PoMoKN* Func_PoMoKN::clone( void ) const
{
    
    return new Func_PoMoKN( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_PoMoKN::createFunction( void ) const
{
    RevBayesCore::TypedDagNode< std::int64_t                          >* na = static_cast<const Natural              &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< std::int64_t                          >* ni = static_cast<const Natural              &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* m  = static_cast<const ModelVector<RealPos> &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* fc = static_cast<const ModelVector<RealPos> &>( this->args[3].getVariable()->getRevObject() ).getDagNode();

    if ( m->getValue().size() != ((na->getValue())*(na->getValue())-na->getValue() ) )
    {
        throw RbException("The number of alleles does not match the number of mutation rates given: n_mut_rates = n_alleles*n_alleles-n_alleles");
    }
    if ( fc->getValue().size() != na->getValue() )
    {
        throw RbException("The number of alleles does not match the number of fitness coefficients given: n_fit_coeff = n_alleles");
    }

    RevBayesCore::PoMoKNRateMatrixFunction* f = new RevBayesCore::PoMoKNRateMatrixFunction( na, ni, m, fc );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_PoMoKN::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "K"       , Natural::getClassTypeSpec(), "Number of alleles", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "N"       , Natural::getClassTypeSpec(), "Number of individuals", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "mu"      , ModelVector<RealPos>::getClassTypeSpec(), "Vector of mutation rates: mu=(mu_a0a1,mu_a1a0,mu_a0a2,mu_a2a0,...)", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "phi"     , ModelVector<RealPos>::getClassTypeSpec(), "Vector of fitness coefficients: phi=(phi_0,phi_1,...,phi_ak)", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );  

        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_PoMoKN::getClassType(void)
{
    
    static std::string rev_type = "Func_PoMoKN";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_PoMoKN::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_PoMoKN::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnPoMoKN";
    
    return f_name;
}


const TypeSpec& Func_PoMoKN::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


