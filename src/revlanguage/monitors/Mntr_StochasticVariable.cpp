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


using namespace RevLanguage;

Mntr_StochasticVariable::Mntr_StochasticVariable(void) : FileMonitor()
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
    const std::string&                  fn      = static_cast<const RlString &> ( filename->getRevObject() ).getValue();
    const std::string&                  sep     = static_cast<const RlString &> ( separator->getRevObject() ).getValue();
    std::uint64_t                       g       = static_cast<const IntegerPos&>( printgen->getRevObject() ).getValue();
    bool                                ap      = static_cast<const RlBoolean &>( append->getRevObject() ).getValue();
    bool                                wv      = static_cast<const RlBoolean &>( version->getRevObject() ).getValue();
    RevBayesCore::StochasticVariableMonitor *m = new RevBayesCore::StochasticVariableMonitor((std::uint64_t)g, fn, sep);
    
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
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        // add the rules from the base class
        const MemberRules &parentRules = FileMonitor::getParameterRules();
        memberRules.insert(memberRules.end(), parentRules.begin(), parentRules.end());
        
        rules_set = true;
    }
    
    return memberRules;
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
    
    FileMonitor::setConstParameter(name, var);
    
}
