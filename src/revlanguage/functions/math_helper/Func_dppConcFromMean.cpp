#include "Func_dppConcFromMean.h"
#include "RbStatisticsHelper.h"
#include "Natural.h"
#include "RealPos.h"

using namespace RevLanguage;

/** default constructor */
Func_dppConcFromMean::Func_dppConcFromMean( void ) : Procedure( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_dppConcFromMean* Func_dppConcFromMean::clone( void ) const
{
    
    return new Func_dppConcFromMean( *this );
}


RevPtr<RevVariable> Func_dppConcFromMean::execute()
{

    double nc = static_cast<const RealPos &>( args[0].getVariable()->getRevObject() ).getValue();
    double ne = static_cast<const Natural &>( args[1].getVariable()->getRevObject() ).getValue();
	
	double meanCP = RevBayesCore::RbStatistics::Helper::dppConcParamFromNumTables(nc, ne);
	RealPos* value = new RealPos( meanCP );

	return new RevVariable( value );
}


/* Get argument rules */
const ArgumentRules& Func_dppConcFromMean::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "numCats"    , RealPos::getClassTypeSpec(), "Number of Categories of the DPP.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "numElements", Natural::getClassTypeSpec(), "Total number of elements.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_dppConcFromMean::getClassType(void)
{
    
    static std::string rev_type = "Func_dppConcFromMean";
    
	return rev_type; 
}

/* Get class type spec describing type of object */
const TypeSpec& Func_dppConcFromMean::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_dppConcFromMean::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnDppConcFromMean";
    
    return f_name;
}


/* Get return type */
const TypeSpec& Func_dppConcFromMean::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = RealPos::getClassTypeSpec();
    
    return return_typeSpec;
}


const TypeSpec& Func_dppConcFromMean::getTypeSpec( void ) const {
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
