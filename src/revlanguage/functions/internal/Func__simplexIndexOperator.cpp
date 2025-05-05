#include "Func__simplexIndexOperator.h"
#include "ModelVector.h"
#include "RlDeterministicNode.h"
#include "RlSimplex.h"
#include "TypedDagNode.h"
#include "VectorIndexOperator.h"

using namespace RevLanguage;

/** default constructor */
Func__simplexIndexOperator::Func__simplexIndexOperator( void ) : TypedFunction<Probability>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func__simplexIndexOperator* Func__simplexIndexOperator::clone( void ) const
{
    
    return new Func__simplexIndexOperator( *this );
}


RevBayesCore::TypedFunction< double >* Func__simplexIndexOperator::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode< RevBayesCore::Simplex >* v = static_cast<const Simplex &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<std::int64_t>* index = static_cast<const Natural &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::VectorIndexOperator<double> *func = new RevBayesCore::VectorIndexOperator<double>(v, index);
    
    return func;
}


/* Get argument rules */
const ArgumentRules& Func__simplexIndexOperator::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "v"    , Simplex::getClassTypeSpec(), "The vector.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "index", Natural::getClassTypeSpec(), "The index.",  ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    
    return argumentRules;
}



const std::string& Func__simplexIndexOperator::getClassType(void)
{
    
    static std::string revClassType = "Func__simplexIndexOperator";
    
    return revClassType;
}

/* Get class type spec describing type of object */
const TypeSpec& Func__simplexIndexOperator::getClassTypeSpec(void)
{
    
    static TypeSpec revClassTypeSpec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return revClassTypeSpec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func__simplexIndexOperator::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "[]";
    
    return f_name;
}



const TypeSpec& Func__simplexIndexOperator::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}

