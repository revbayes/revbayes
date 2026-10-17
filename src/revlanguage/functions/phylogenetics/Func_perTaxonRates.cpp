#include "Func_perTaxonRates.h"

#include "RevVariable.h"
#include "RlDeterministicNode.h"
#include "TypeSpec.h"

using namespace RevLanguage;


Func_perTaxonRates::Func_perTaxonRates( void ) : Func_FossilRateEstimator()
{

}


Func_perTaxonRates* Func_perTaxonRates::clone( void ) const
{
    return new Func_perTaxonRates( *this );
}


const std::string& Func_perTaxonRates::getClassType( void )
{
    static std::string rev_type = "Func_perTaxonRates";

    return rev_type;
}


const TypeSpec& Func_perTaxonRates::getClassTypeSpec( void )
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), &Func_FossilRateEstimator::getClassTypeSpec() );

    return rev_type_spec;
}


/** Get the Rev name of the function */
std::string Func_perTaxonRates::getFunctionName( void ) const
{
    static std::string f_name = "fnPerTaxonRates";

    return f_name;
}


const TypeSpec& Func_perTaxonRates::getTypeSpec( void ) const
{
    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}


RevBayesCore::FossilRateEstimatorFunction::Method Func_perTaxonRates::method( void ) const
{
    return RevBayesCore::FossilRateEstimatorFunction::PER_TAXON;
}
