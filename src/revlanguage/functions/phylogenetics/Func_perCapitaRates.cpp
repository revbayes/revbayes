#include "Func_perCapitaRates.h"

#include "RevVariable.h"
#include "RlDeterministicNode.h"
#include "TypeSpec.h"

using namespace RevLanguage;


Func_perCapitaRates::Func_perCapitaRates( void ) : Func_FossilRateEstimator()
{

}


Func_perCapitaRates* Func_perCapitaRates::clone( void ) const
{
    return new Func_perCapitaRates( *this );
}


const std::string& Func_perCapitaRates::getClassType( void )
{
    static std::string rev_type = "Func_perCapitaRates";

    return rev_type;
}


const TypeSpec& Func_perCapitaRates::getClassTypeSpec( void )
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), &Func_FossilRateEstimator::getClassTypeSpec() );

    return rev_type_spec;
}


/** Get the Rev name of the function */
std::string Func_perCapitaRates::getFunctionName( void ) const
{
    static std::string f_name = "fnPerCapitaRates";

    return f_name;
}


const TypeSpec& Func_perCapitaRates::getTypeSpec( void ) const
{
    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}


RevBayesCore::FossilRateEstimatorFunction::Method Func_perCapitaRates::method( void ) const
{
    return RevBayesCore::FossilRateEstimatorFunction::PER_TAXON;
}
