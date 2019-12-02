#include "Func_FreeBinary.h"
#include "FreeBinaryRateMatrixFunction.h"
#include "FreeBinaryRateMatrixFunction.h"
#include "ModelVector.h"
#include "Natural.h"
#include "RateMatrix_JC.h"
#include "Real.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "RlRateMatrix.h"
#include "RlSimplex.h"
#include "RlBoolean.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_FreeBinary::Func_FreeBinary( void ) : TypedFunction<RateMatrix>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_FreeBinary* Func_FreeBinary::clone( void ) const
{
    
    return new Func_FreeBinary( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_FreeBinary::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* glr = static_cast<const ModelVector<RealPos> &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    bool r = static_cast<const RlBoolean &>( this->args[1].getVariable()->getRevObject() ).getDagNode()->getValue();

    RevBayesCore::FreeBinaryRateMatrixFunction* f = new RevBayesCore::FreeBinaryRateMatrixFunction( glr, r );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_FreeBinary::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "transition_rates", ModelVector<Real>::getClassTypeSpec(), "The transition rates between the two states.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "rescaled",          RlBoolean::getClassTypeSpec(),              "Should the matrix be normalized?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_FreeBinary::getClassType(void)
{
    
    static std::string rev_type = "Func_FreeBinary";
    
	return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_FreeBinary::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_FreeBinary::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnFreeBinary";
    
    return f_name;
}


const TypeSpec& Func_FreeBinary::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
