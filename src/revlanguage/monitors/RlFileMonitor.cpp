#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RlBoolean.h"
#include "RlFileMonitor.h"
#include "IntegerPos.h"
#include "RlString.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "TypeSpec.h"

namespace RevBayesCore { class DagNode; }

using namespace RevLanguage;

FileMonitor::FileMonitor(void) : Monitor()
{
    
}


const std::string& FileMonitor::getClassType(void)
{
    
    static std::string rev_type = "FileMonitor";
    
	return rev_type; 
}


const TypeSpec& FileMonitor::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Monitor::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/** Return member rules (no members) */
const MemberRules& FileMonitor::getParameterRules(void) const
{

    static MemberRules memberRules;
    static bool rules_set = false;

    if ( !rules_set )
    {
        memberRules.push_back( new ArgumentRule("append"    , RlBoolean::getClassTypeSpec() , "Should we append or overwrite if the file exists?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );
        memberRules.push_back( new ArgumentRule("filename"  , RlString::getClassTypeSpec()  , "The name of the file for storing the samples.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule("printgen"  , IntegerPos::getClassTypeSpec(), "The number of generations between stored samples.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new IntegerPos(1) ) );
        std::vector<std::string> sep = {"separator", "delimiter"};
        memberRules.push_back( new ArgumentRule(sep , RlString::getClassTypeSpec()          , "The separator/delimiter between columns in the file.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString("\t") ) );
        memberRules.push_back( new ArgumentRule("version"   , RlBoolean::getClassTypeSpec() , "Should we record the software version?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );

        rules_set = true;
    }

    return memberRules;
}


/** Set a member variable */
void FileMonitor::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{

    if ( name == "filename" )
    {
        filename = var;
    }
    else if ( name == "separator" || name == "delimiter" || name == "separator/delimiter" )
    {
        separator = var;
    }
    else if ( name == "printgen" )
    {
        printgen = var;
    }
    else if ( name == "append" )
    {
        append = var;
    }
    else if ( name == "version" )
    {
        version = var;
    }
    else
    {
        Monitor::setConstParameter(name, var);
    }

}




