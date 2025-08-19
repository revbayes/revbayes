#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RlMonitor.h"
#include "Argument.h"
#include "MemberProcedure.h"
#include "MethodTable.h"
#include "Monitor.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlUtils.h"
#include "StringUtilities.h"
#include "TypeSpec.h"
#include "WorkspaceToCoreWrapperObject.h"

namespace RevBayesCore { class DagNode; }

using namespace RevLanguage;

Monitor::Monitor(void) : WorkspaceToCoreWrapperObject<RevBayesCore::Monitor>()
{
    
    // add method for call "addVariable" as a function
    ArgumentRules* addArgRules = new ArgumentRules();
    addArgRules->push_back( new ArgumentRule( "x" , RevObject::getClassTypeSpec(), "A variable you want to monitor.", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberProcedure( "addVariable", RlUtils::Void, addArgRules) );
    
}


Monitor::Monitor(RevBayesCore::Monitor *m) : WorkspaceToCoreWrapperObject<RevBayesCore::Monitor>( m )
{
    
    // add method for call "addVariable" as a function
    ArgumentRules* addArgRules = new ArgumentRules();
    addArgRules->push_back( new ArgumentRule( "x" , RevObject::getClassTypeSpec(), "A variable you want to monitor.", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberProcedure( "addVariable", RlUtils::Void, addArgRules) );
    
}



RevPtr<RevVariable> Monitor::executeMethod(const std::string& name, const std::vector<Argument>& args, bool &found)
{
    
    if ( name == "addVariable" )
    {
        found = true;
        
        RevBayesCore::DagNode* node = args[0].getVariable()->getRevObject().getDagNode();
        
        RevBayesCore::Monitor *m = static_cast<RevBayesCore::Monitor*>(this->value);
        m->addVariable( node );
        
        return NULL;
    }
    else if ( name == "removeVariable" )
    {
        found = true;
        
        RevBayesCore::DagNode* node = args[0].getVariable()->getRevObject().getDagNode();
        
        RevBayesCore::Monitor *m = static_cast<RevBayesCore::Monitor*>(this->value);
        m->removeVariable( node );
        
        return NULL;
    }
    
    return WorkspaceToCoreWrapperObject<RevBayesCore::Monitor>::executeMethod( name, args, found );
}


const std::string& Monitor::getClassType(void)
{
    
    static std::string rev_type = "Monitor";
    
	return rev_type; 
}


const TypeSpec& Monitor::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( WorkspaceToCoreWrapperObject<RevBayesCore::Monitor>::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the aliases of the Rev name for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Monitor::getConstructorFunctionAliases( void ) const
{
    // create a constructor function name alias variable that is the same for all instance of this class
    std::vector<std::string> aliases;
    
    std::vector<std::string> monitor_aliases = getMonitorAliases();
    for (size_t i=0; i < monitor_aliases.size(); ++i)
    {
        std::string tmp = monitor_aliases[i];
        std::string c_name = "mn" + StringUtilities::firstCharToUpper( tmp );
        
        aliases.push_back( c_name );
    }
    
    return aliases;
}


/**
 * Get the Rev name for the constructor function (with appropriate prefix).
 *
 * \return Rev name of constructor function.
 */
std::string Monitor::getConstructorFunctionName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string tmp = getMonitorName();
    std::string c_name = "mn" + StringUtilities::firstCharToUpper( tmp );
    
    return c_name;
}


void Monitor::printValue(std::ostream &o, bool user) const
{    
    printValue(o);
}


void Monitor::printValue(std::ostream &o) const
{
    o << getMonitorName();
}




