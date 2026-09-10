#include "Func_PoMo4N.h"
#include "PoMo4NRateMatrixFunction.h"
#include "RateMatrix_PoMo4N.h"
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
Func_PoMo4N::Func_PoMo4N( void ) : TypedFunction<RateMatrix>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_PoMo4N* Func_PoMo4N::clone( void ) const
{
    
    return new Func_PoMo4N( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_PoMo4N::createFunction( void ) const
{
    RevBayesCore::TypedDagNode< std::int64_t                          >* ni = static_cast<const Natural              &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* m  = static_cast<const ModelVector<RealPos> &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* fc = static_cast<const ModelVector<RealPos> &>( this->args[2].getVariable()->getRevObject() ).getDagNode();

    if ( m->getValue().size() != 12 )
    {
        throw RbException("The number of mutation rates should match the number of pairwise combinations of alleles: 12.");
    }
    if ( fc->getValue().size() != 4 )
    {
        throw RbException("The number of fitness coefficients given should match the number of alleles: 4.");
    }

    RevBayesCore::PoMo4NRateMatrixFunction* f = new RevBayesCore::PoMo4NRateMatrixFunction( ni, m, fc );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_PoMo4N::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "N"       , Natural::getClassTypeSpec(), "Number of individuals", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "mu"      , ModelVector<RealPos>::getClassTypeSpec(), "Vector of mutation rates: mu=(mu_a0a1,mu_a1a0,mu_a0a2,mu_a2a0...)", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "phi"     , ModelVector<RealPos>::getClassTypeSpec(), "Vector of fitness coefficients: phi=(phi_0,phi_1,phi_2,phi_3)", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );  

        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_PoMo4N::getClassType(void)
{
    
    static std::string rev_type = "Func_PoMo4N";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_PoMo4N::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_PoMo4N::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnPoMo4N";
    
    return f_name;
}


const TypeSpec& Func_PoMo4N::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


