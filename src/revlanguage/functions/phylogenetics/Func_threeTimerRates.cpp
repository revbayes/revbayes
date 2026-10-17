#include "Func_threeTimerRates.h"

#include "RevVariable.h"
#include "RlDeterministicNode.h"
#include "TypeSpec.h"

using namespace RevLanguage;


Func_threeTimerRates::Func_threeTimerRates( void ) : Func_FossilRateEstimator()
{

}


Func_threeTimerRates* Func_threeTimerRates::clone( void ) const
{
    return new Func_threeTimerRates( *this );
}


const std::string& Func_threeTimerRates::getClassType( void )
{
    static std::string rev_type = "Func_threeTimerRates";

    return rev_type;
}


const TypeSpec& Func_threeTimerRates::getClassTypeSpec( void )
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), &Func_FossilRateEstimator::getClassTypeSpec() );

    return rev_type_spec;
}


/** Get the Rev name of the function */
std::string Func_threeTimerRates::getFunctionName( void ) const
{
    static std::string f_name = "fnThreeTimerRates";

    return f_name;
}


const TypeSpec& Func_threeTimerRates::getTypeSpec( void ) const
{
    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}


RevBayesCore::FossilRateEstimatorFunction::Method Func_threeTimerRates::method( void ) const
{
    return RevBayesCore::FossilRateEstimatorFunction::THREE_TIMER;
}
