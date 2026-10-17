#include "Func_boundaryCrosserRates.h"

#include "RevVariable.h"
#include "RlDeterministicNode.h"
#include "TypeSpec.h"

using namespace RevLanguage;


Func_boundaryCrosserRates::Func_boundaryCrosserRates( void ) : Func_FossilRateEstimator()
{

}


Func_boundaryCrosserRates* Func_boundaryCrosserRates::clone( void ) const
{
    return new Func_boundaryCrosserRates( *this );
}


const std::string& Func_boundaryCrosserRates::getClassType( void )
{
    static std::string rev_type = "Func_boundaryCrosserRates";

    return rev_type;
}


const TypeSpec& Func_boundaryCrosserRates::getClassTypeSpec( void )
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), &Func_FossilRateEstimator::getClassTypeSpec() );

    return rev_type_spec;
}


/** Get the Rev name of the function */
std::string Func_boundaryCrosserRates::getFunctionName( void ) const
{
    static std::string f_name = "fnBoundaryCrosserRates";

    return f_name;
}


const TypeSpec& Func_boundaryCrosserRates::getTypeSpec( void ) const
{
    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}


RevBayesCore::FossilRateEstimatorFunction::Method Func_boundaryCrosserRates::method( void ) const
{
    return RevBayesCore::FossilRateEstimatorFunction::BOUNDARY_CROSSER;
}
