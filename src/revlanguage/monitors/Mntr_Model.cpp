#include <ostream>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Mntr_Model.h"
#include "ModelMonitor.h"
#include "ModelVector.h"
#include "Natural.h"
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

Mntr_Model::Mntr_Model(void) : FileMonitor()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Mntr_Model* Mntr_Model::clone(void) const 
{
    
	return new Mntr_Model(*this);
}


void Mntr_Model::constructInternalObject( void ) 
{
    // we free the memory first
    delete value;
    
    // now allocate a new sliding move
    const std::string&                  fn      = static_cast<const RlString &>( filename->getRevObject() ).getValue();
    const std::string&                  sep     = static_cast<const RlString &>( separator->getRevObject() ).getValue();
    unsigned int                        g       = (int)static_cast<const IntegerPos  &>( printgen->getRevObject() ).getValue();
    bool                                pp      = static_cast<const RlBoolean &>( posterior->getRevObject() ).getValue();
    bool                                l       = static_cast<const RlBoolean &>( likelihood->getRevObject() ).getValue();
    bool                                pr      = static_cast<const RlBoolean &>( prior->getRevObject() ).getValue();
    bool                                ap      = static_cast<const RlBoolean &>( append->getRevObject() ).getValue();
    bool                                so      = static_cast<const RlBoolean &>( stochOnly->getRevObject() ).getValue();
    bool                                wv      = static_cast<const RlBoolean &>( version->getRevObject() ).getValue();
    const std::string&                 fmt      = static_cast<const RlString &>( format->getRevObject() ).getValue();

    ModelVector<RlString> excl = static_cast<const ModelVector<RlString> &>(exclude->getRevObject());
    std::set<std::string> exclude_list;
    for ( size_t i = 0; i < excl.size(); ++i ) {
        exclude_list.insert(excl[i]);
    }

    SampleFormat Format = (fmt == "json") ? SampleFormat(JSONFormat()) : SampleFormat(SeparatorFormat(sep));
    RevBayesCore::ModelMonitor *m = new RevBayesCore::ModelMonitor((unsigned long)g, fn, Format, exclude_list);
    
    // now set the flags
    m->setAppend( ap );
    m->setPrintLikelihood( l );
    m->setPrintPosterior( pp );
    m->setPrintPrior( pr );
    m->setPrintVersion( wv );
    m->setStochasticNodesOnly( so );
    
    // store the new model into our value variable
    value = m;
}


/** Get Rev type of object */
const std::string& Mntr_Model::getClassType(void) 
{ 
    
    static std::string rev_type = "Mntr_Model";
    
	return rev_type; 
}

/** Get class type spec describing type of object */
const TypeSpec& Mntr_Model::getClassTypeSpec(void) 
{ 
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Monitor::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Mntr_Model::getMonitorName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "Model";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Mntr_Model::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set ) 
    {
        
        memberRules.push_back( new ArgumentRule("posterior"     , RlBoolean::getClassTypeSpec(), "Should we print the joint posterior probability?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );
        memberRules.push_back( new ArgumentRule("likelihood"    , RlBoolean::getClassTypeSpec(), "Should we print the likelihood?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );
        memberRules.push_back( new ArgumentRule("prior"         , RlBoolean::getClassTypeSpec(), "Should we print the joint prior probability?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );
        memberRules.push_back( new ArgumentRule("stochasticOnly", RlBoolean::getClassTypeSpec(), "Should we monitor stochastic variables only?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );
        memberRules.push_back( new ArgumentRule{"exclude", ModelVector<RlString>::getClassTypeSpec(), "Variables to exclude from the monitor", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new ModelVector<RlString>()});
        memberRules.push_back( new ArgumentRule("format"        , RlString::getClassTypeSpec(),  "Output format", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString("separator") ) );
        
        // add the rules from the base class
        const MemberRules &parentRules = FileMonitor::getParameterRules();
        memberRules.insert(memberRules.end(), parentRules.begin(), parentRules.end());

        rules_set = true;
    }
    
    return memberRules;
}

/** Get type spec */
const TypeSpec& Mntr_Model::getTypeSpec( void ) const 
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void Mntr_Model::printValue(std::ostream &o) const 
{
    
    o << "Mntr_Model";
}


/** Set a member variable */
void Mntr_Model::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) 
{
    
    if ( name == "prior" )
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
    else if ( name == "stochasticOnly" ) 
    {
        stochOnly = var;
    }
    else if ( name == "exclude" )
    {
        exclude = var;
    }
    else if ( name == "format" )
    {
        format = var;
    }
    else 
    {
        FileMonitor::setConstParameter(name, var);
    }
    
}
