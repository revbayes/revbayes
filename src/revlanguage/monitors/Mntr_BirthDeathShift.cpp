#include "Mntr_BirthDeathShift.h"

#include <cstddef>
#include <string>

#include "ArgumentRule.h"
#include "IntegerPos.h"
#include "RevObject.h"
#include "RlTimeTree.h"
#include "RlString.h"
#include "FastBirthDeathShiftProcess.h"
#include "TypeSpec.h"
#include "ArgumentRules.h"
#include "Monitor.h"
#include "RbBoolean.h"
#include "RlBoolean.h"
#include "RlTree.h"
#include "BirthDeathShiftMonitor.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore { class Tree; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }


using namespace RevLanguage;


Mntr_BirthDeathShift::Mntr_BirthDeathShift(void) : FileMonitor()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Mntr_BirthDeathShift* Mntr_BirthDeathShift::clone(void) const
{
    
    return new Mntr_BirthDeathShift(*this);
}


void Mntr_BirthDeathShift::constructInternalObject( void )
{
    
    const std::string& file_name      = static_cast<const RlString  &>( filename->getRevObject()           ).getValue();
    const std::string& sep            = static_cast<const RlString  &>( separator->getRevObject()          ).getValue();
    unsigned int       print_gen      = (int)static_cast<const IntegerPos   &>( printgen->getRevObject()      ).getValue();
    bool               app            = static_cast<const RlBoolean &>( append->getRevObject()             ).getValue();
    bool               wv             = static_cast<const RlBoolean &>( version->getRevObject()            ).getValue();
    
    RevBayesCore::BirthDeathShiftMonitor *m;

    if ( static_cast<const RevLanguage::Tree&>( bdsp->getRevObject() ).isModelObject() )
    {
        RevBayesCore::TypedDagNode<RevBayesCore::Tree>* bdsp_tdn = static_cast<const RevLanguage::Tree&>( bdsp->getRevObject() ).getDagNode();
        RevBayesCore::StochasticNode<RevBayesCore::Tree>* bdsp_sn  = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree>* >( bdsp_tdn );

        RevBayesCore::FastBirthDeathShiftProcess *sse_process = NULL;
        sse_process = dynamic_cast<RevBayesCore::FastBirthDeathShiftProcess*>( &bdsp_sn->getDistribution() );
        sse_process->setSampleCharacterHistory( true );

        m = new RevBayesCore::BirthDeathShiftMonitor( bdsp_sn, (unsigned long)print_gen, file_name, sep );
        m->setAppend( app );
        m->setPrintVersion( wv );
    }
    else
    {
    	throw RbException("Must provide a FastBirthDeathShift object.");
    }
    
    delete value;
    value = m;
}


/** Get Rev type of object */
const std::string& Mntr_BirthDeathShift::getClassType(void)
{
    
    static std::string revType = "Mntr_BirthDeathShift";
    
    return revType;
    
}


/** Get class type spec describing type of object */
const TypeSpec& Mntr_BirthDeathShift::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Monitor::getClassTypeSpec() ) );
    
    return rev_type_spec;
    
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Mntr_BirthDeathShift::getMonitorName( void ) const
{
    
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "BirthDeathShift";
    
    return c_name;
    
}


/** Return member rules (no members) */
const MemberRules& Mntr_BirthDeathShift::getParameterRules(void) const
{
    
    static MemberRules monitor_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        monitor_rules.push_back( new ArgumentRule("bdsp"          , TimeTree::getClassTypeSpec(),  "The birth-death-shift process to monitor.",                      ArgumentRule::BY_REFERENCE, ArgumentRule::ANY, NULL ) );
        // add the rules from the base class
        const MemberRules &parentRules = FileMonitor::getParameterRules();
        monitor_rules.insert(monitor_rules.end(), parentRules.begin(), parentRules.end());

        rules_set = true;
    }
    
    return monitor_rules;
    
}


/** Get type spec */
const TypeSpec& Mntr_BirthDeathShift::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
    
}


/** Get type spec */
void Mntr_BirthDeathShift::printValue(std::ostream &o) const
{
    
    o << "Mntr_BirthDeathShift";
    
}


/** Set a member variable */
void Mntr_BirthDeathShift::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "bdsp" )
    {
        bdsp = var;
    }
    else
    {
        FileMonitor::setConstParameter(name, var);
    }
    
}
