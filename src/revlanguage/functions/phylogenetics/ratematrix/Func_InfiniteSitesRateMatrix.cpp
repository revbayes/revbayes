#include "Func_InfiniteSitesRateMatrix.h"

#include <cstddef>

#include "InfiniteSitesRateMatrixFunction.h"
#include "Natural.h"
#include "RlRateMatrix.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "TypeSpec.h"

using namespace RevLanguage;

/** default constructor */
Func_InfiniteSitesRateMatrix::Func_InfiniteSitesRateMatrix( void ) : TypedFunction<RateMatrix>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_InfiniteSitesRateMatrix* Func_InfiniteSitesRateMatrix::clone( void ) const
{
    
    return new Func_InfiniteSitesRateMatrix( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_InfiniteSitesRateMatrix::createFunction( void ) const
{
    
    int ns = (int)static_cast<const Natural &>( this->args[0].getVariable()->getRevObject() ).getValue();
    RevBayesCore::InfiniteSitesRateMatrixFunction* f = new RevBayesCore::InfiniteSitesRateMatrixFunction( size_t(ns) );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_InfiniteSitesRateMatrix::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "num_states", Natural::getClassTypeSpec(), "The number of states.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(2) ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_InfiniteSitesRateMatrix::getClassType(void)
{
    
    static std::string rev_type = "Func_InfiniteSitesRateMatrix";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_InfiniteSitesRateMatrix::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_InfiniteSitesRateMatrix::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnInfiniteSites";
    
    return f_name;
}


const TypeSpec& Func_InfiniteSitesRateMatrix::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
