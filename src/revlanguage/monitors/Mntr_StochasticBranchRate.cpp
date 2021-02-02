#include "Mntr_StochasticBranchRate.h"

#include <stddef.h>
#include <string>

#include "ArgumentRule.h"
#include "IntegerPos.h"
#include "RlMonitor.h"
#include "RevObject.h"
#include "RlTimeTree.h"
#include "RlString.h"
#include "StateDependentSpeciationExtinctionProcess.h"
#include "TypeSpec.h"
#include "ArgumentRules.h"
#include "Monitor.h"
#include "RbBoolean.h"
#include "RlBoolean.h"
#include "RlTree.h"
#include "StochasticBranchRateMonitor.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore { class Tree; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }


using namespace RevLanguage;


Mntr_StochasticBranchRate::Mntr_StochasticBranchRate(void) : Monitor()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Mntr_StochasticBranchRate* Mntr_StochasticBranchRate::clone(void) const
{
    
    return new Mntr_StochasticBranchRate(*this);
}


void Mntr_StochasticBranchRate::constructInternalObject( void )
{
    
    const std::string& file_name      = static_cast<const RlString  &>( filename->getRevObject()           ).getValue();
    const std::string& sep            = static_cast<const RlString  &>( separator->getRevObject()          ).getValue();
    unsigned int                print_gen      = (int)static_cast<const IntegerPos   &>( printgen->getRevObject()      ).getValue();
    bool               app            = static_cast<const RlBoolean &>( append->getRevObject()             ).getValue();
    bool               wv             = static_cast<const RlBoolean &>( version->getRevObject()            ).getValue();
    
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* cdbdp_tdn = static_cast<const RevLanguage::Tree&>( cdbdp->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::Tree>* cdbdp_sn  = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree>* >( cdbdp_tdn );
    
    RevBayesCore::StateDependentSpeciationExtinctionProcess *sse_process = NULL;
    sse_process = dynamic_cast<RevBayesCore::StateDependentSpeciationExtinctionProcess*>( &cdbdp_sn->getDistribution() );
    sse_process->setSampleCharacterHistory( true );
    
    RevBayesCore::StochasticBranchRateMonitor *m = new RevBayesCore::StochasticBranchRateMonitor( cdbdp_sn, (unsigned long)print_gen, file_name, sep );
    m->setAppend( app );
    m->setPrintVersion( wv );
    
    delete value;
    value = m;
    
}


/** Get Rev type of object */
const std::string& Mntr_StochasticBranchRate::getClassType(void)
{
    
    static std::string revType = "Mntr_StochasticBranchRate";
    
    return revType;
    
}


/** Get class type spec describing type of object */
const TypeSpec& Mntr_StochasticBranchRate::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Monitor::getClassTypeSpec() ) );
    
    return rev_type_spec;
    
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Mntr_StochasticBranchRate::getMonitorName( void ) const
{
    
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "StochasticBranchRate";
    
    return c_name;
    
}


/** Return member rules (no members) */
const MemberRules& Mntr_StochasticBranchRate::getParameterRules(void) const
{
    
    static MemberRules monitor_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        monitor_rules.push_back( new ArgumentRule("cdbdp"          , TimeTree::getClassTypeSpec(),  "The character dependent birth-death process to monitor.",                      ArgumentRule::BY_REFERENCE, ArgumentRule::ANY ) );
        monitor_rules.push_back( new ArgumentRule("filename"       , RlString::getClassTypeSpec() , "The file to save sampled character histories.",                                ArgumentRule::BY_VALUE,     ArgumentRule::ANY ) );
        monitor_rules.push_back( new ArgumentRule("printgen"       , IntegerPos::getClassTypeSpec()  , "How frequently (in number of iterations) should we save sampled character histories? 1 by default.",                              ArgumentRule::BY_VALUE,     ArgumentRule::ANY, new IntegerPos(1) ) );
        monitor_rules.push_back( new ArgumentRule("separator"      , RlString::getClassTypeSpec() , "The delimiter between variables. \t by default.",                              ArgumentRule::BY_VALUE,     ArgumentRule::ANY, new RlString("\t") ) );
        monitor_rules.push_back( new ArgumentRule("append"         , RlBoolean::getClassTypeSpec(), "Should we append to an existing file? False by default.",                  ArgumentRule::BY_VALUE,     ArgumentRule::ANY, new RlBoolean(false) ) );
        monitor_rules.push_back( new ArgumentRule("version"        , RlBoolean::getClassTypeSpec(), "Should we record the software version?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );
        
        rules_set = true;
    }
    
    return monitor_rules;
    
}


/** Get type spec */
const TypeSpec& Mntr_StochasticBranchRate::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
    
}


/** Get type spec */
void Mntr_StochasticBranchRate::printValue(std::ostream &o) const
{
    
    o << "Mntr_StochasticBranchRate";
    
}


/** Set a member variable */
void Mntr_StochasticBranchRate::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "" )
    {
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
    else if ( name == "cdbdp" )
    {
        cdbdp = var;
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
