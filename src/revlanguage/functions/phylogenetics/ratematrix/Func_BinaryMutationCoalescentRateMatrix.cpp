#include "BinaryMutationCoalescentRateMatrixFunction.h"
#include "Func_BinaryMutationCoalescentRateMatrix.h"
#include "ModelVector.h"
#include "Natural.h"
#include "RateMatrix_BinaryMutationCoalescent.h"
#include "Real.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "RlRateMatrix.h"
#include "RlSimplex.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_BinaryMutationCoalescentRateMatrix::Func_BinaryMutationCoalescentRateMatrix( void ) : TypedFunction<RateMatrix>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_BinaryMutationCoalescentRateMatrix* Func_BinaryMutationCoalescentRateMatrix::clone( void ) const
{
    
    return new Func_BinaryMutationCoalescentRateMatrix( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_BinaryMutationCoalescentRateMatrix::createFunction( void ) const
{
    
    size_t n = static_cast<const Natural &>( this->args[0].getVariable()->getRevObject() ).getValue();

    RevBayesCore::TypedDagNode< double >* mu = static_cast<const RealPos &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< double >* ne = static_cast<const RealPos &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    
    RevBayesCore::BinaryMutationCoalescentRateMatrixFunction* f = new RevBayesCore::BinaryMutationCoalescentRateMatrixFunction( n, mu, ne );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_BinaryMutationCoalescentRateMatrix::getArgumentRules( void ) const
{
    
    static ArgumentRules argument_rules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( rules_set == false )
    {
        
        argument_rules.push_back( new ArgumentRule( "n",  Natural::getClassTypeSpec(), "The number of samples.",         ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "mu", RealPos::getClassTypeSpec(), "The mutation rate.",             ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "Ne", RealPos::getClassTypeSpec(), "The effective population size.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argument_rules;
}


const std::string& Func_BinaryMutationCoalescentRateMatrix::getClassType(void)
{
    
    static std::string rev_type = "Func_BinaryMutationCoalescentRateMatrix";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_BinaryMutationCoalescentRateMatrix::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_BinaryMutationCoalescentRateMatrix::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnBinaryMutationCoalescentRateMatrix";
    
    return f_name;
}


const TypeSpec& Func_BinaryMutationCoalescentRateMatrix::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
