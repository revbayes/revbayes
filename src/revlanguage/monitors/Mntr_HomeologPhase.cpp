#include "Mntr_HomeologPhase.h"

#include <stddef.h>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RlMonitor.h"
#include "IntegerPos.h"
#include "RevObject.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlString.h"
#include "TypeSpec.h"
#include "HomeologPhaseMonitor.h"
#include "Monitor.h"
#include "RbBoolean.h"
#include "RlBoolean.h"
#include "StochasticNode.h"

namespace RevBayesCore { class AbstractHomologousDiscreteCharacterData; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }

using namespace RevLanguage;

Mntr_HomeologPhase::Mntr_HomeologPhase(void) : Monitor()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Mntr_HomeologPhase* Mntr_HomeologPhase::clone(void) const
{
    
    return new Mntr_HomeologPhase(*this);
}


void Mntr_HomeologPhase::constructInternalObject( void )
{
    const std::string&                  fn      = static_cast<const RlString &>( filename->getRevObject() ).getValue();
    const std::string&                  sep     = static_cast<const RlString &>( separator->getRevObject() ).getValue();
    unsigned int                                 g       = (int)static_cast<const IntegerPos  &>( printgen->getRevObject() ).getValue();
    
    bool                                ap      = static_cast<const RlBoolean &>( append->getRevObject() ).getValue();
    bool                                wv      = static_cast<const RlBoolean &>( version->getRevObject() ).getValue();
    

    
    RevBayesCore::TypedDagNode<RevBayesCore::AbstractHomologousDiscreteCharacterData>* ctmc_tdn = NULL;
    RevBayesCore::StochasticNode<RevBayesCore::AbstractHomologousDiscreteCharacterData>* ctmc_sn = NULL;
    
    ctmc_tdn = static_cast<const RevLanguage::AbstractHomologousDiscreteCharacterData&>( ctmc->getRevObject() ).getDagNode();
    ctmc_sn  = static_cast<RevBayesCore::StochasticNode<RevBayesCore::AbstractHomologousDiscreteCharacterData>* >(ctmc_tdn);
        
    
    
    delete value;
    RevBayesCore::HomeologPhaseMonitor *m;
    m = new RevBayesCore::HomeologPhaseMonitor(ctmc_sn, (unsigned long)g, fn, sep);
    m->setAppend( ap );
    m->setPrintVersion(wv);
    value = m;
    
}


/** Get Rev type of object */
const std::string& Mntr_HomeologPhase::getClassType(void)
{
    
    static std::string rev_type = "Mntr_HomeologPhase";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Mntr_HomeologPhase::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Monitor::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Mntr_HomeologPhase::getMonitorName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "HomeologPhase";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Mntr_HomeologPhase::getParameterRules(void) const
{
    
    static MemberRules asMonitorMemberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        asMonitorMemberRules.push_back( new ArgumentRule("ctmc"           , AbstractHomologousDiscreteCharacterData::getClassTypeSpec(), "", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY, NULL ) );
        asMonitorMemberRules.push_back( new ArgumentRule("filename"       , RlString::getClassTypeSpec() , "", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        asMonitorMemberRules.push_back( new ArgumentRule("printgen"       , IntegerPos::getClassTypeSpec()  , "", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new IntegerPos(1) ) );
        asMonitorMemberRules.push_back( new ArgumentRule("separator"      , RlString::getClassTypeSpec() , "", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString("\t") ) );
        asMonitorMemberRules.push_back( new ArgumentRule("append"         , RlBoolean::getClassTypeSpec(), "", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );
        asMonitorMemberRules.push_back( new ArgumentRule("version"        , RlBoolean::getClassTypeSpec(), "Should we record the software version?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );

        rules_set = true;
    }
    
    return asMonitorMemberRules;
}

/** Get type spec */
const TypeSpec& Mntr_HomeologPhase::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void Mntr_HomeologPhase::printValue(std::ostream &o) const
{
    
    o << "Mntr_HomeologPhase";
}


/** Set a member variable */
void Mntr_HomeologPhase::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "" ) {
        vars.push_back( var );
    }
    else if ( name == "filename" )
    {
        filename = var;
    }
    else if ( name == "separator" )
    {
        separator = var;
    }
    else if ( name == "ctmc" )
    {
        ctmc = var;
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
