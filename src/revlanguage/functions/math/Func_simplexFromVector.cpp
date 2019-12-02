#include "ArgumentRule.h"
#include "Func_simplexFromVector.h"
#include "ModelVector.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "RlSimplex.h"
#include "SimplexFromVectorFunction.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"


using namespace RevLanguage;


/** Default constructor */
Func_simplexFromVector::Func_simplexFromVector( void ) :
    TypedFunction<Simplex>()
{
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_simplexFromVector* Func_simplexFromVector::clone( void ) const
{
    return new Func_simplexFromVector( *this );
}


/** Execute function: Compute simplex from vector of RealPos values */
RevBayesCore::TypedFunction< RevBayesCore::Simplex >* Func_simplexFromVector::createFunction( void ) const
{
    
    const RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* vec;
    vec = static_cast< const ModelVector<RealPos>& >( args[0].getVariable()->getRevObject() ).getDagNode();
    
    RevBayesCore::SimplexFromVectorFunction*    func    = new RevBayesCore::SimplexFromVectorFunction( vec );
    
    return func;
}


/** Get argument rules */
const ArgumentRules& Func_simplexFromVector::getArgumentRules( void ) const
{
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "x", ModelVector<RealPos>::getClassTypeSpec(), "A vector of numbers.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object (static version) */
const std::string& Func_simplexFromVector::getClassType( void )
{
    static std::string rev_type = "Func_simplexFromVector";
    
	return rev_type;
}


/** Get Rev type spec of object (static version) */
const TypeSpec& Func_simplexFromVector::getClassTypeSpec( void )
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), &Function::getClassTypeSpec() );
    
	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_simplexFromVector::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "simplex";
    
    return f_name;
}


/** Get Rev type spec from an instance of the object */
const TypeSpec& Func_simplexFromVector::getTypeSpec( void ) const
{
    return getClassTypeSpec();
}

