//
//  Mntr_StochasticCharacterMapping.cpp
//  RevBayes_development_branch
//
//  Created by will freyman on 12/13/16.
//  Copyright (c) 2016 will freyman. All rights reserved.
//

#include "Mntr_StochasticCharacterMapping.h"

#include "AbstractCharacterData.h"
#include "ArgumentRule.h"
#include "Model.h"
#include "Natural.h"
#include "NaturalNumbersState.h"
#include "RlMonitor.h"
#include "RbException.h"
#include "RevObject.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlModel.h"
#include "RlTimeTree.h"
#include "RlBranchLengthTree.h"
#include "RlString.h"
#include "StandardState.h"
#include "StateDependentSpeciationExtinctionProcess.h"
#include "TypeSpec.h"


using namespace RevLanguage;


Mntr_StochasticCharacterMapping::Mntr_StochasticCharacterMapping(void) : Monitor()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Mntr_StochasticCharacterMapping* Mntr_StochasticCharacterMapping::clone(void) const
{
    
    return new Mntr_StochasticCharacterMapping(*this);
}


void Mntr_StochasticCharacterMapping::constructInternalObject( void )
{
    
    const std::string& file_name      = static_cast<const RlString  &>( filename->getRevObject()           ).getValue();
    bool               is             = static_cast<const RlBoolean &>( include_simmap->getRevObject()     ).getValue();
    bool               sd             = static_cast<const RlBoolean &>( use_simmap_default->getRevObject() ).getValue();
    const std::string& sep            = static_cast<const RlString  &>( separator->getRevObject()          ).getValue();
    int                print_gen      = (int)static_cast<const Natural   &>( printgen->getRevObject()      ).getValue();
    bool               app            = static_cast<const RlBoolean &>( append->getRevObject()             ).getValue();
    bool               wv             = static_cast<const RlBoolean &>( version->getRevObject()            ).getValue();
    
    
    RevBayesCore::TypedDagNode<RevBayesCore::AbstractHomologousDiscreteCharacterData>* ctmc_tdn = NULL;
    RevBayesCore::StochasticNode<RevBayesCore::AbstractHomologousDiscreteCharacterData>* ctmc_sn = NULL;
    
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* cdbdp_tdn = NULL;
    RevBayesCore::StochasticNode<RevBayesCore::Tree>* cdbdp_sn = NULL;
    
    if ( static_cast<const RevLanguage::AbstractHomologousDiscreteCharacterData&>( ctmc->getRevObject() ).isModelObject() )
    {
        ctmc_tdn = static_cast<const RevLanguage::AbstractHomologousDiscreteCharacterData&>( ctmc->getRevObject() ).getDagNode();
        ctmc_sn  = static_cast<RevBayesCore::StochasticNode<RevBayesCore::AbstractHomologousDiscreteCharacterData>* >(ctmc_tdn);
    }
    else if ( static_cast<const RevLanguage::Tree&>( cdbdp->getRevObject() ).isModelObject() )
    {
        cdbdp_tdn = static_cast<const RevLanguage::Tree&>( cdbdp->getRevObject() ).getDagNode();
        cdbdp_sn  = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree>* >( cdbdp_tdn );
        
        RevBayesCore::StateDependentSpeciationExtinctionProcess *sse_process = NULL;
        sse_process = dynamic_cast<RevBayesCore::StateDependentSpeciationExtinctionProcess*>( &cdbdp_sn->getDistribution() );
        sse_process->setSampleCharacterHistory( true );
    }
    else
    {
        throw RbException("mnStochasticCharacterMap requires either a CTMC or a character-dependent birth death process (CDBDP).");
    }
    
    
    RevBayesCore::StochasticCharacterMappingMonitor<RevBayesCore::NaturalNumbersState> *m;
    if ( static_cast<const RevLanguage::AbstractHomologousDiscreteCharacterData&>( ctmc->getRevObject() ).isModelObject() )
    {
        m = new RevBayesCore::StochasticCharacterMappingMonitor<RevBayesCore::NaturalNumbersState>( ctmc_sn, (unsigned long)print_gen, file_name, is, sd, sep );
    }
    else
    {
        m = new RevBayesCore::StochasticCharacterMappingMonitor<RevBayesCore::NaturalNumbersState>( cdbdp_sn, (unsigned long)print_gen, file_name, is, sd, sep );
    }
    m->setAppend( app );
    m->setPrintVersion( wv );
    
    delete value;
    value = m;
    
}


/** Get Rev type of object */
const std::string& Mntr_StochasticCharacterMapping::getClassType(void)
{
    
    static std::string revType = "Mntr_StochasticCharacterMap";
    
    return revType;
    
}


/** Get class type spec describing type of object */
const TypeSpec& Mntr_StochasticCharacterMapping::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Monitor::getClassTypeSpec() ) );
    
    return rev_type_spec;
    
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Mntr_StochasticCharacterMapping::getMonitorName( void ) const
{
    
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "StochasticCharacterMap";
    
    return c_name;
    
}


/** Return member rules (no members) */
const MemberRules& Mntr_StochasticCharacterMapping::getParameterRules(void) const
{
    
    static MemberRules monitor_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        monitor_rules.push_back( new ArgumentRule("ctmc"           , AbstractHomologousDiscreteCharacterData::getClassTypeSpec(), "The continuous-time Markov process to monitor.", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY, NULL ) );
        monitor_rules.push_back( new ArgumentRule("cdbdp"          , TimeTree::getClassTypeSpec(),  "The character dependent birth-death process to monitor.",                      ArgumentRule::BY_REFERENCE, ArgumentRule::ANY, NULL) );
        monitor_rules.push_back( new ArgumentRule("filename"       , RlString::getClassTypeSpec() , "The file to save sampled character histories.",                                ArgumentRule::BY_VALUE,     ArgumentRule::ANY ) );
        monitor_rules.push_back( new ArgumentRule("include_simmap" , RlBoolean::getClassTypeSpec(), "Should we log SIMMAP/phytools compatible newick strings? True by default.",    ArgumentRule::BY_VALUE,     ArgumentRule::ANY, new RlBoolean(true) ) );
        monitor_rules.push_back( new ArgumentRule("use_simmap_default" , RlBoolean::getClassTypeSpec(), "Should we use the default SIMMAP/phytools event ordering? True by default.",    ArgumentRule::BY_VALUE,     ArgumentRule::ANY, new RlBoolean(true) ) );
        monitor_rules.push_back( new ArgumentRule("printgen"       , Natural::getClassTypeSpec()  , "How frequently (in number of iterations) should we save sampled character histories? 1 by default.",                              ArgumentRule::BY_VALUE,     ArgumentRule::ANY, new Natural(1) ) );
        monitor_rules.push_back( new ArgumentRule("separator"      , RlString::getClassTypeSpec() , "The delimiter between variables. \t by default.",                              ArgumentRule::BY_VALUE,     ArgumentRule::ANY, new RlString("\t") ) );
        monitor_rules.push_back( new ArgumentRule("append"         , RlBoolean::getClassTypeSpec(), "Should we append to an existing file? False by default.",                  ArgumentRule::BY_VALUE,     ArgumentRule::ANY, new RlBoolean(false) ) );
        monitor_rules.push_back( new ArgumentRule("version"        , RlBoolean::getClassTypeSpec(), "Should we record the software version?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );

        rules_set = true;
    }
    
    return monitor_rules;
    
}


/** Get type spec */
const TypeSpec& Mntr_StochasticCharacterMapping::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
    
}


/** Get type spec */
void Mntr_StochasticCharacterMapping::printValue(std::ostream &o) const
{
    
    o << "Mntr_StochasticCharacterMap";
    
}


/** Set a member variable */
void Mntr_StochasticCharacterMapping::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
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
    else if ( name == "cdbdp" )
    {
        cdbdp = var;
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
    else if ( name == "include_simmap" )
    {
        include_simmap = var;
    }
    else if ( name == "use_simmap_default" )
    {
        use_simmap_default = var;
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
