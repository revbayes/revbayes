#include "Real.h"
#include "RlString.h"
#include "RlUtils.h"

using namespace RevLanguage;


RevObject* RlUtils::RlTypeConverter::toReal(double x) 
{
    return new Real( x );
}

RevObject* RlUtils::RlTypeConverter::toString(const std::string &x) 
{
    return new RlString( x );
}
