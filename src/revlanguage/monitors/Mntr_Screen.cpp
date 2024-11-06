
#include <algorithm>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Ellipsis.h"
#include "Mntr_Screen.h"
#include "Natural.h"
#include "IntegerPos.h"
#include "ScreenMonitor.h"
#include "TypeSpec.h"
#include "Monitor.h"
#include "RbBoolean.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlBoolean.h"
#include "RlMonitor.h"

namespace RevBayesCore { class DagNode; }


using namespace RevLanguage;

Mntr_Screen::Mntr_Screen(void) : Monitor()
{}


Mntr_Screen* Mntr_Screen::clone(void) const
{    
	return new Mntr_Screen(*this);
}


void Mntr_Screen::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate space for a new Mntr_Screen object
    std::uint64_t g = (unsigned long) static_cast<const IntegerPos &>( printgen->getRevObject() ).getValue();
    vars.erase( unique( vars.begin(), vars.end() ), vars.end() );
    sort( vars.begin(), vars.end(), compareVarNames );
    std::vector<RevBayesCore::DagNode *> n;
    for (std::vector<RevPtr<const RevVariable> >::iterator i = vars.begin(); i != vars.end(); ++i)
    {
        RevBayesCore::DagNode* node = (*i)->getRevObject().getDagNode();
        n.push_back( node );
    }
    bool pp = static_cast<const RlBoolean &>( posterior->getRevObject() ).getValue();
    bool l = static_cast<const RlBoolean &>( likelihood->getRevObject() ).getValue();
    bool pr = static_cast<const RlBoolean &>( prior->getRevObject() ).getValue();
    value = new RevBayesCore::ScreenMonitor(n, g, pp, l, pr);
}


const std::string& Mntr_Screen::getClassType(void)
{   
    static std::string rev_type = "Mntr_Screen";    
	return rev_type; 
}


const TypeSpec& Mntr_Screen::getClassTypeSpec(void)
{    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Monitor::getClassTypeSpec() ) );   
	return rev_type_spec; 
}


std::string Mntr_Screen::getMonitorName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "Screen";    
    return c_name;
}


const MemberRules& Mntr_Screen::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        memberRules.push_back( new Ellipsis( "Variables to monitor.", RevObject::getClassTypeSpec() ) );
        memberRules.push_back( new ArgumentRule("printgen"  , IntegerPos::getClassTypeSpec(), "The number of generations between displayed samples.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new IntegerPos(1) ) );
        memberRules.push_back( new ArgumentRule("posterior" , RlBoolean::getClassTypeSpec() , "Monitor the joint posterior probability.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );
        memberRules.push_back( new ArgumentRule("likelihood", RlBoolean::getClassTypeSpec() , "Monitor the joint likelihood.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );
        memberRules.push_back( new ArgumentRule("prior"     , RlBoolean::getClassTypeSpec() , "Monitor the joint prior probability.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );
        
        rules_set = true;
    }
    
    return memberRules;
}


const TypeSpec& Mntr_Screen::getTypeSpec( void ) const
{    
    static TypeSpec type_spec = getClassTypeSpec();    
    return type_spec;
}


void Mntr_Screen::printValue(std::ostream &o) const
{  
    o << "Mntr_Screen";
}


void Mntr_Screen::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "" )
    {
        vars.push_back( var );
    }
    else if ( name == "printgen" )
    {   
        printgen = var;
    }
    else if ( name == "prior" )
    {
        prior = var;
    }
    else if ( name == "posterior" )
    {
        posterior = var;
    }
    else if ( name == "likelihood" )
    {
        likelihood = var;
    }
    else
    {
        RevObject::setConstParameter(name, var);
    }
}
