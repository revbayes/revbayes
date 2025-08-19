
#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "DemographicFunction.h"
#include "RlDemographicFunction.h"
#include "TypeSpec.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "StringUtilities.h"
#include "WorkspaceToCoreWrapperObject.h"

namespace RevLanguage { class Argument; }


using namespace RevLanguage;

DemographicFunction::DemographicFunction(void) : WorkspaceToCoreWrapperObject<RevBayesCore::DemographicFunction>()
{
    
}


DemographicFunction::DemographicFunction(RevBayesCore::DemographicFunction *m) : WorkspaceToCoreWrapperObject<RevBayesCore::DemographicFunction>( m )
{
    
}



RevPtr<RevVariable> DemographicFunction::executeMethod(const std::string& name, const std::vector<Argument>& args, bool &found)
{
    
    return WorkspaceToCoreWrapperObject<RevBayesCore::DemographicFunction>::executeMethod( name, args, found );
}


/** Get Rev type of object */
const std::string& DemographicFunction::getClassType(void)
{
    
    static std::string rev_type = "DemographicFunction";
    
	return rev_type; 
}

/** Get class type spec describing type of object */
const TypeSpec& DemographicFunction::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( WorkspaceToCoreWrapperObject<RevBayesCore::DemographicFunction>::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the aliases of the Rev name for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> DemographicFunction::getConstructorFunctionAliases( void ) const
{
    // create a constructor function name alias variable that is the same for all instance of this class
    std::vector<std::string> aliases;
    
    std::vector<std::string> DemographicFunction_aliases = getDemographicFunctionAliases();
    for (size_t i=0; i < DemographicFunction_aliases.size(); ++i)
    {
        std::string tmp = DemographicFunction_aliases[i];
        std::string c_name = "df" + StringUtilities::firstCharToUpper( tmp );
        
        aliases.push_back( c_name );
    }
    
    return aliases;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string DemographicFunction::getConstructorFunctionName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string tmp = getDemographicFunctionName();
    std::string c_name = "df" + StringUtilities::firstCharToUpper( tmp );
    
    return c_name;
}


/** Get type spec */
void DemographicFunction::printValue(std::ostream &o, bool user) const
{
    
    printValue(o);
}

/** Get type spec */
void DemographicFunction::printValue(std::ostream &o) const
{
    
    o << getDemographicFunctionName();
}




