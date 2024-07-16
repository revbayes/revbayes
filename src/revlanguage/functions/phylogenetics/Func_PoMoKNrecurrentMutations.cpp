#include "Func_PoMoKNrecurrentMutations.h"
#include "PoMoKNrecurrentMutationsRateMatrixFunction.h"
#include "RateMatrix_PoMoKNrecurrentMutations.h"
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
Func_PoMoKNrecurrentMutations::Func_PoMoKNrecurrentMutations( void ) : TypedFunction<RateMatrix>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_PoMoKNrecurrentMutations* Func_PoMoKNrecurrentMutations::clone( void ) const
{
    
    return new Func_PoMoKNrecurrentMutations( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_PoMoKNrecurrentMutations::createFunction( void ) const
{
    long n_alleles  = static_cast<const Natural                                                                                     &>( this->args[0].getVariable()->getRevObject() ).getValue();
    long n_virtual  = static_cast<const Natural                                                                                     &>( this->args[1].getVariable()->getRevObject() ).getValue();

    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* mutation_rates          = static_cast<const ModelVector<RealPos>   &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* fitness_coefficients    = static_cast<const ModelVector<RealPos>   &>( this->args[3].getVariable()->getRevObject() ).getDagNode();

    bool recurrent_mutations = static_cast<const RlBoolean                                                                          &>( this->args[4].getVariable()->getRevObject() ).getValue();



    if ( mutation_rates->getValue().size() != (n_alleles*n_alleles-n_alleles ) )
    {
        throw RbException("The number of alleles does not match the number of mutation rates given: n_mut_rates = n_alleles*n_alleles-n_alleles");
    }
    if ( fitness_coefficients->getValue().size() != n_alleles )
    {
        throw RbException("The number of alleles does not match the number of fitness coefficients given: n_fit_coeff = n_alleles");
    }

    RevBayesCore::PoMoKNrecurrentMutationsRateMatrixFunction* f = new RevBayesCore::PoMoKNrecurrentMutationsRateMatrixFunction( n_alleles, n_virtual , mutation_rates, fitness_coefficients, recurrent_mutations );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_PoMoKNrecurrentMutations::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "K"    , Natural::getClassTypeSpec(), "Number of alleles", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "V"    , Natural::getClassTypeSpec(), "Number of virtual individuals", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "mu"   , ModelVector<RealPos>::getClassTypeSpec(), "Vector of mutation rates: mu=(mu_a0a1,mu_a1a0,mu_a0a2,mu_a2a0,...)", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "phi"  , ModelVector<RealPos>::getClassTypeSpec(), "Vector of fitness coefficients: phi=(phi_0,phi_1,...,phi_ak)", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );  
        argumentRules.push_back( new ArgumentRule( "rec"  , RlBoolean::getClassTypeSpec(), "Mutation model: TRUE for recurrent and FALSE for boundary mutations", ArgumentRule::BY_VALUE, ArgumentRule::ANY,  new RlBoolean( true )  ) );

        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_PoMoKNrecurrentMutations::getClassType(void)
{
    
    static std::string rev_type = "Func_PoMoKNrecurrentMutations";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_PoMoKNrecurrentMutations::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_PoMoKNrecurrentMutations::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnPoMoKNmutations";
    
    return f_name;
}


const TypeSpec& Func_PoMoKNrecurrentMutations::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


