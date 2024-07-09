#include "Func_revPoMo2Nrecurrent.h"
#include "revPoMo2NrecurrentRateMatrixFunction.h"
#include "RateMatrix_revPoMo2Nrecurrent.h"
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
Func_revPoMo2Nrecurrent::Func_revPoMo2Nrecurrent( void ) : TypedFunction<RateMatrix>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_revPoMo2Nrecurrent* Func_revPoMo2Nrecurrent::clone( void ) const
{
    
    return new Func_revPoMo2Nrecurrent( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_revPoMo2Nrecurrent::createFunction( void ) const
{
    RevBayesCore::TypedDagNode< long                          >* ni = static_cast<const Natural              &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::Simplex          >* bf = static_cast<const Simplex              &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< double                        >* ex = static_cast<const RealPos              &>( this->args[2].getVariable()->getRevObject() ).getDagNode();

    if ( bf->getValue().size() != 2 )
    {
        throw RbException("This model only admits two alleles and therefore only two allele frequencies are expected.");
    }

    RevBayesCore::revPoMo2NrecurrentRateMatrixFunction* f = new RevBayesCore::revPoMo2NrecurrentRateMatrixFunction( ni, bf, ex );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_revPoMo2Nrecurrent::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "N"   , Natural::getClassTypeSpec(), "Number of individuals", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "pi"  , Simplex::getClassTypeSpec(), "Vector of allele frequencies: pi=(pi_a0,pi_a1)", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "rho" , RealPos::getClassTypeSpec(), "Exchangeability: rho_a0a1", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_revPoMo2Nrecurrent::getClassType(void)
{
    
    static std::string rev_type = "Func_revPoMo2Nrecurrent";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_revPoMo2Nrecurrent::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_revPoMo2Nrecurrent::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnRecurrentMutationsPoMo2N";
    
    return f_name;
}


const TypeSpec& Func_revPoMo2Nrecurrent::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


