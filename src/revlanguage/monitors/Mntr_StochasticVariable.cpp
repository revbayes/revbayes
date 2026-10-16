#include "RevVariableUtilities.h"
#include <ostream>
#include "TypeUtilities.h"
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
    const std::string&                  fn      = valueFor<RlString>(filename);
    const std::string&                  sep     = valueFor<RlString>(separator);
    std::uint64_t                       g       = valueFor<IntegerPos>(printgen);
    bool                                ap      = valueFor<RlBoolean>(append);
    bool                                wv      = valueFor<RlBoolean>(version);
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
