//
//  Mntr_BirthDeathStochasticCharacterMapping.cpp
//  RevBayes_development_branch
//
//  Created by will freyman on 12/13/16.
//  Copyright (c) 2016 will freyman. All rights reserved.
//

#include "Mntr_BirthDeathStochasticCharacterMapping.h"

#include <stddef.h>
#include <string>

#include "ArgumentRule.h"
#include "BirthDeathStochasticCharacterMappingMonitor.h"
#include "Natural.h"
#include "NaturalNumbersState.h"
#include "StandardState.h"
#include "RlMonitor.h"
#include "RbException.h"
#include "RevObject.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlTimeTree.h"
#include "RlString.h"
//#include "StateDependentSpeciationExtinctionProcess.h"
#include "TypeSpec.h"
#include "ArgumentRules.h"
#include "DiscreteTaxonData.h"
#include "Monitor.h"
#include "RbBoolean.h"
#include "RlBoolean.h"
#include "RlTree.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore { class AbstractHomologousDiscreteCharacterData; }
namespace RevBayesCore { class Tree; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }


using namespace RevLanguage;


Mntr_BirthDeathStochasticCharacterMapping::Mntr_BirthDeathStochasticCharacterMapping(void) : Monitor()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Mntr_BirthDeathStochasticCharacterMapping* Mntr_BirthDeathStochasticCharacterMapping::clone(void) const
{
    
    return new Mntr_BirthDeathStochasticCharacterMapping(*this);
}


void Mntr_BirthDeathStochasticCharacterMapping::constructInternalObject( void )
{
    
    const std::string& file_name      = static_cast<const RlString        &>( filename->getRevObject()           ).getValue();
    bool               is             = static_cast<const RlBoolean       &>( include_simmap->getRevObject()     ).getValue();
    bool               sd             = static_cast<const RlBoolean       &>( use_simmap_default->getRevObject() ).getValue();
    const std::string& sep            = static_cast<const RlString        &>( separator->getRevObject()          ).getValue();
    int                print_gen      = (int)static_cast<const Natural    &>( printgen->getRevObject()           ).getValue();
    bool               app            = static_cast<const RlBoolean       &>( append->getRevObject()             ).getValue();
    bool               wv             = static_cast<const RlBoolean       &>( version->getRevObject()            ).getValue();
    size_t             idx            = (size_t)static_cast<const Natural &>( index->getRevObject()              ).getValue();
    
    
    RevBayesCore::TypedDagNode<RevBayesCore::AbstractHomologousDiscreteCharacterData>* ctmc_tdn = NULL;
    RevBayesCore::StochasticNode<RevBayesCore::AbstractHomologousDiscreteCharacterData>* ctmc_sn = NULL;
    
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tree_tdn = NULL;
    RevBayesCore::StochasticNode<RevBayesCore::Tree>* tree_sn = NULL;
    
//    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* cdbdp_tdn = NULL;
//    RevBayesCore::StochasticNode<RevBayesCore::Tree>* cdbdp_sn = NULL;
    
    if ( static_cast<const RevLanguage::AbstractHomologousDiscreteCharacterData&>( ctmc->getRevObject() ).isModelObject() )
    {
        ctmc_tdn = static_cast<const RevLanguage::AbstractHomologousDiscreteCharacterData&>( ctmc->getRevObject() ).getDagNode();
        ctmc_sn  = static_cast<RevBayesCore::StochasticNode<RevBayesCore::AbstractHomologousDiscreteCharacterData>* >(ctmc_tdn);
        
        tree_tdn = static_cast<const RevLanguage::TimeTree&>( tree->getRevObject() ).getDagNode();
        tree_sn  = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree>* >(tree_tdn);
    }
//    else if ( static_cast<const RevLanguage::Tree&>( cdbdp->getRevObject() ).isModelObject() )
//    {
//        
//        throw RbException("mnStochasticCharacterMap does not currently support `complete_tree=true` character-dependent birth death process (CDBDP) models.");
//        
//        cdbdp_tdn = static_cast<const RevLanguage::Tree&>( cdbdp->getRevObject() ).getDagNode();
//        cdbdp_sn  = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree>* >( cdbdp_tdn );
//        
//        RevBayesCore::StateDependentSpeciationExtinctionProcess *sse_process = NULL;
//        sse_process = dynamic_cast<RevBayesCore::StateDependentSpeciationExtinctionProcess*>( &cdbdp_sn->getDistribution() );
//        sse_process->setSampleCharacterHistory( true );
//    }
    else
    {
        throw RbException("mnStochasticCharacterMap requires either a CTMC or a character-dependent birth death process (CDBDP).");
    }
    
    
    std::string model_type = "";
    std::string data_type = "";
    if ( static_cast<const RevLanguage::AbstractHomologousDiscreteCharacterData&>( ctmc->getRevObject() ).isModelObject() ) {
        model_type = "ctmc";
        data_type = ctmc_sn->getValue().getDataType();
    }
    else {
        model_type = "cdbdp";
        data_type = "NaturalNumbers";
    }
    
    
    if ( model_type == "ctmc" )
    {
        std::string data_type = ctmc_sn->getValue().getDataType();
        if (data_type == "Standard") {
            
            RevBayesCore::BirthDeathStochasticCharacterMappingMonitor<RevBayesCore::StandardState>* m;
            m = new RevBayesCore::BirthDeathStochasticCharacterMappingMonitor<RevBayesCore::StandardState>( ctmc_sn, tree_sn, (unsigned long)print_gen, file_name, is, sd, sep, idx - 1 );
            m->setAppend( app );
            m->setPrintVersion( wv );
            
            delete value;
            value = m;
            
        } else if (data_type == "NaturalNumbers") {
            RevBayesCore::BirthDeathStochasticCharacterMappingMonitor<RevBayesCore::NaturalNumbersState>* m;
            m = new RevBayesCore::BirthDeathStochasticCharacterMappingMonitor<RevBayesCore::NaturalNumbersState>( ctmc_sn, tree_sn, (unsigned long)print_gen, file_name, is, sd, sep, idx - 1 );
            m->setAppend( app );
            m->setPrintVersion( wv );
            
            delete value;
            value = m;
        }
        
    }
//    else if ( model_type == "cdbdp" )
//    {
//        if (data_type == "NaturalNumbers") {
//            RevBayesCore::BirthDeathStochasticCharacterMappingMonitor<RevBayesCore::NaturalNumbersState>* m;
//            m = new RevBayesCore::BirthDeathStochasticCharacterMappingMonitor<RevBayesCore::NaturalNumbersState>( cdbdp_sn, (unsigned long)print_gen, file_name, is, sd, sep );
//            m->setAppend( app );
//            m->setPrintVersion( wv );
//            
//            delete value;
//            value = m;
//        }
//    }
    
    
}


/** Get Rev type of object */
const std::string& Mntr_BirthDeathStochasticCharacterMapping::getClassType(void)
{
    
    static std::string revType = "Mntr_StochasticCharacterMap";
    
    return revType;
    
}


/** Get class type spec describing type of object */
const TypeSpec& Mntr_BirthDeathStochasticCharacterMapping::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Monitor::getClassTypeSpec() ) );
    
    return rev_type_spec;
    
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Mntr_BirthDeathStochasticCharacterMapping::getMonitorName( void ) const
{
    
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "BirthDeathStochasticCharacterMap";
    
    return c_name;
    
}


/** Return member rules (no members) */
const MemberRules& Mntr_BirthDeathStochasticCharacterMapping::getParameterRules(void) const
{
    
    static MemberRules monitor_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        monitor_rules.push_back( new ArgumentRule("ctmc"           , AbstractHomologousDiscreteCharacterData::getClassTypeSpec(), "The continuous-time Markov process to monitor.", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY) );
        monitor_rules.push_back( new ArgumentRule("tree"           , TimeTree::getClassTypeSpec(),  "The birth-death process to generate subclades.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC) );
        monitor_rules.push_back( new ArgumentRule("filename"       , RlString::getClassTypeSpec() , "The file to save sampled character histories.",                                ArgumentRule::BY_VALUE,     ArgumentRule::ANY ) );
        monitor_rules.push_back( new ArgumentRule("include_simmap" , RlBoolean::getClassTypeSpec(), "Should we log SIMMAP/phytools compatible newick strings? True by default.",    ArgumentRule::BY_VALUE,     ArgumentRule::ANY, new RlBoolean(true) ) );
        monitor_rules.push_back( new ArgumentRule("use_simmap_default" , RlBoolean::getClassTypeSpec(), "Should we use the default SIMMAP/phytools event ordering? True by default.",    ArgumentRule::BY_VALUE,     ArgumentRule::ANY, new RlBoolean(true) ) );
        monitor_rules.push_back( new ArgumentRule("printgen"       , Natural::getClassTypeSpec()  , "How frequently (in number of iterations) should we save sampled character histories? 1 by default.",                              ArgumentRule::BY_VALUE,     ArgumentRule::ANY, new Natural(1) ) );
        monitor_rules.push_back( new ArgumentRule("separator"      , RlString::getClassTypeSpec() , "The delimiter between variables. \t by default.",                              ArgumentRule::BY_VALUE,     ArgumentRule::ANY, new RlString("\t") ) );
        monitor_rules.push_back( new ArgumentRule("append"         , RlBoolean::getClassTypeSpec(), "Should we append to an existing file? False by default.",                  ArgumentRule::BY_VALUE,     ArgumentRule::ANY, new RlBoolean(false) ) );
        monitor_rules.push_back( new ArgumentRule("version"        , RlBoolean::getClassTypeSpec(), "Should we record the software version?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );
        monitor_rules.push_back( new ArgumentRule("index"          , Natural::getClassTypeSpec(), "The index of the character to be mointored.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(1) ) );
        
        rules_set = true;
    }
    
    return monitor_rules;
    
}


/** Get type spec */
const TypeSpec& Mntr_BirthDeathStochasticCharacterMapping::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
    
}


/** Get type spec */
void Mntr_BirthDeathStochasticCharacterMapping::printValue(std::ostream &o) const
{
    
    o << "Mntr_BirthDeathStochasticCharacterMap";
    
}


/** Set a member variable */
void Mntr_BirthDeathStochasticCharacterMapping::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
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
    else if ( name == "tree" )
    {
        tree = var;
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
    else if ( name == "index" )
    {
        index = var;
    }
    else
    {
        Monitor::setConstParameter(name, var);
    }
    
}
