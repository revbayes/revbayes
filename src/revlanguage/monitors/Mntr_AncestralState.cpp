#include "Mntr_AncestralState.h"

#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "AncestralStateMonitor.h"
#include "RlMonitor.h"
#include "Natural.h"
#include "RbException.h"
#include "RevObject.h"
#include "RlString.h"
#include "TypeSpec.h"
#include "NaturalNumbersState.h"
#include "DnaState.h"
#include "Monitor.h"
#include "RbBoolean.h"
#include "RlBoolean.h"
#include "RlTree.h"
#include "StandardState.h"
#include "TypedDistribution.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class Tree; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }

using namespace RevLanguage;

Mntr_AncestralState::Mntr_AncestralState(void) : Monitor()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Mntr_AncestralState* Mntr_AncestralState::clone(void) const
{
    
    return new Mntr_AncestralState(*this);
}


void Mntr_AncestralState::constructInternalObject( void )
{
    const std::string&                  fn      = static_cast<const RlString &>( filename->getRevObject() ).getValue();
    const std::string&                  sep     = static_cast<const RlString &>( separator->getRevObject() ).getValue();
    int                                 g       = (int)static_cast<const Natural  &>( printgen->getRevObject() ).getValue();
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* t = static_cast<const Tree &>( tree->getRevObject() ).getDagNode();
    RevBayesCore::DagNode*				ch		= ctmc->getRevObject().getDagNode();
    bool                                ap      = static_cast<const RlBoolean &>( append->getRevObject() ).getValue();
    bool                                wv      = static_cast<const RlBoolean &>( version->getRevObject() ).getValue();
    std::string							character = static_cast<const RlString &>( monitorType->getRevObject() ).getValue();
    
    delete value;
    if (character == "NaturalNumbers")
    {
        
        RevBayesCore::AncestralStateMonitor<RevBayesCore::NaturalNumbersState> *m = new RevBayesCore::AncestralStateMonitor<RevBayesCore::NaturalNumbersState>(t, ch, (unsigned long)g, fn, sep);
        m->setAppend( ap );
        m->setPrintVersion( wv );
        value = m;
        
    }
    else if (character == "DNA")
    {
        
        RevBayesCore::AncestralStateMonitor<RevBayesCore::DnaState> *m = new RevBayesCore::AncestralStateMonitor<RevBayesCore::DnaState>(t, ch, (unsigned long)g, fn, sep);
        m->setAppend( ap );
        m->setPrintVersion( wv );
        value = m;
        
    }
    else if (character == "StandardState")
    {
        
        RevBayesCore::AncestralStateMonitor<RevBayesCore::StandardState> *m = new RevBayesCore::AncestralStateMonitor<RevBayesCore::StandardState>(t, ch, (unsigned long)g, fn, sep);
        m->setAppend( ap );
        m->setPrintVersion( wv );
        value = m;
        
    }
    else
    {
        throw RbException( "Incorrect character type specified. Valid options are: NaturalNumbers, DNA" );
    }
    
    
    
    
}


/** Get Rev type of object */
const std::string& Mntr_AncestralState::getClassType(void)
{
    
    static std::string rev_type = "Mntr_AncestralState";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& Mntr_AncestralState::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Monitor::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Mntr_AncestralState::getMonitorName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "AncestralState";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Mntr_AncestralState::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        memberRules.push_back( new ArgumentRule("tree"          , Tree::getClassTypeSpec()     , "The tree which we monitor.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule("ctmc"          , RevObject::getClassTypeSpec(), "The CTMC process.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        memberRules.push_back( new ArgumentRule("filename"      , RlString::getClassTypeSpec() , "The name of the file for storing the samples.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule("type"          , RlString::getClassTypeSpec() , "The type of data to store.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule("printgen"      , Natural::getClassTypeSpec()  , "The frequency how often to sample.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(1) ) );
        memberRules.push_back( new ArgumentRule("separator"     , RlString::getClassTypeSpec() , "The separator between columns in the file.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString("\t") ) );
        memberRules.push_back( new ArgumentRule("append"        , RlBoolean::getClassTypeSpec(), "Should we append or overwrite if the file exists?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );
        memberRules.push_back( new ArgumentRule("version"   , RlBoolean::getClassTypeSpec(), "Should we record the software version?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );

        rules_set = true;
    }
    
    return memberRules;
}

/** Get type spec */
const TypeSpec& Mntr_AncestralState::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void Mntr_AncestralState::printValue(std::ostream &o) const
{
    
    o << "Mntr_AncestralState";
}


/** Set a member variable */
void Mntr_AncestralState::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
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
    else if ( name == "tree" ) {
        tree = var;
    }
    else if ( name == "type" )
    {
        monitorType = var;
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
