#include <ostream>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Mntr_StochasticVariable.h"
#include "StochasticVariableMonitor.h"
#include "IntegerPos.h"
#include "RevObject.h"
#include "RlString.h"
#include "TypeSpec.h"
#include "Monitor.h"
#include "RbBoolean.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlBoolean.h"
#include "RlMonitor.h"


using namespace RevLanguage;

Mntr_StochasticVariable::Mntr_StochasticVariable(void) : Monitor()
{

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Mntr_StochasticVariable* Mntr_StochasticVariable::clone(void) const
{
    
    return new Mntr_StochasticVariable(*this);
}


void Mntr_StochasticVariable::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new sliding move
    const std::string&                  fn      = static_cast<const RlString &>( filename->getRevObject() ).getValue();
    const std::string&                  sep     = static_cast<const RlString &>( separator->getRevObject() ).getValue();
    unsigned long                                g       = static_cast<const IntegerPos  &>( printgen->getRevObject() ).getValue();
    bool                                ap      = static_cast<const RlBoolean &>( append->getRevObject() ).getValue();
    bool                                wv      = static_cast<const RlBoolean &>( version->getRevObject() ).getValue();
    RevBayesCore::StochasticVariableMonitor *m = new RevBayesCore::StochasticVariableMonitor((unsigned long)g, fn, sep);
    
    // now set the flags
    m->setAppend( ap );
    m->setPrintVersion( wv );
    
    // store the new StochasticVariable into our value variable
    value = m;
}


/** Get Rev type of object */
const std::string& Mntr_StochasticVariable::getClassType(void)
{
    
    static std::string rev_type = "Mntr_StochasticVariable";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& Mntr_StochasticVariable::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Monitor::getClassTypeSpec() ) );
    
    return rev_type_spec;
}



/** Return member rules (no members) */
const MemberRules& Mntr_StochasticVariable::getParameterRules(void) const
{
    
    static MemberRules StochasticVariableMonitorMemberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        StochasticVariableMonitorMemberRules.push_back( new ArgumentRule("filename"  , RlString::getClassTypeSpec() , "The name of the file.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        StochasticVariableMonitorMemberRules.push_back( new ArgumentRule("printgen"  , IntegerPos::getClassTypeSpec()  , "The frequency how often we print.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new IntegerPos(1) ) );
        StochasticVariableMonitorMemberRules.push_back( new ArgumentRule("separator" , RlString::getClassTypeSpec() , "The delimiter between variables.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString("\t") ) );
        StochasticVariableMonitorMemberRules.push_back( new ArgumentRule("append"    , RlBoolean::getClassTypeSpec(), "Should we append or overwrite if the file exists?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );
        StochasticVariableMonitorMemberRules.push_back( new ArgumentRule("version"   , RlBoolean::getClassTypeSpec(), "Should we record the software version?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );

        
        rules_set = true;
    }
    
    return StochasticVariableMonitorMemberRules;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Mntr_StochasticVariable::getMonitorName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "StochasticVariable";
    
    return c_name;
}


/** Get type spec */
const TypeSpec& Mntr_StochasticVariable::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void Mntr_StochasticVariable::printValue(std::ostream &o) const
{
    
    o << "Mntr_StochasticVariable";
}


/** Set a member variable */
void Mntr_StochasticVariable::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "filename" )
    {
        filename = var;
    }
    else if ( name == "separator" )
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
