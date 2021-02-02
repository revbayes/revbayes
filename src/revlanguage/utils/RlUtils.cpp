#include <iosfwd>

#include "Real.h"
#include "RlString.h"
#include "RlUtils.h"

namespace RevLanguage { class RevObject; }

using namespace RevLanguage;


RevObject* RlUtils::RlTypeConverter::toReal(double x) 
{
    return new Real( x );
}

RevObject* RlUtils::RlTypeConverter::toString(const std::string &x) 
{
    return new RlString( x );
}
