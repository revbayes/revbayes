
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Ellipsis.h"
#include "ModelVector.h"
#include "Mntr_Screen.h"
#include "Natural.h"
#include "RbException.h"
#include "RevObject.h"
#include "RlString.h"
#include "ScreenMonitor.h"
#include "TypeSpec.h"


using namespace RevLanguage;

Mntr_Screen::Mntr_Screen(void) : Monitor()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Mntr_Screen* Mntr_Screen::clone(void) const
{
    
	return new Mntr_Screen(*this);
}


void Mntr_Screen::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate space for a new Mntr_Screen object
    int g = (int)static_cast<const Natural &>( printgen->getRevObject() ).getValue();

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


/** Get Rev type of object */
const std::string& Mntr_Screen::getClassType(void)
{
    
    static std::string rev_type = "Mntr_Screen";
    
	return rev_type; 
}

/** Get class type spec describing type of object */
const TypeSpec& Mntr_Screen::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Monitor::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Mntr_Screen::getMonitorName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "Screen";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Mntr_Screen::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        memberRules.push_back( new Ellipsis( "Variables to monitor.", RevObject::getClassTypeSpec() ) );
        memberRules.push_back( new ArgumentRule("printgen"  , Natural::getClassTypeSpec()  , "The frequency how often the variables are monitored.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(1) ) );
        memberRules.push_back( new ArgumentRule("posterior" , RlBoolean::getClassTypeSpec(), "Monitor the joint posterior probability.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );
        memberRules.push_back( new ArgumentRule("likelihood", RlBoolean::getClassTypeSpec(), "Monitor the joint likelihood.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );
        memberRules.push_back( new ArgumentRule("prior"     , RlBoolean::getClassTypeSpec(), "Monitor the joint prior probability.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );
        
        rules_set = true;
    }
    
    return memberRules;
}

/** Get type spec */
const TypeSpec& Mntr_Screen::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void Mntr_Screen::printValue(std::ostream &o) const
{
    
    o << "Mntr_Screen";
}


/** Set a member variable */
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
