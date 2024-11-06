//
//  Mntr_StochasticBranchStateTimes.cpp
//  RevBayes_development_branch
//
//  Created by will freyman on 12/13/16.
//  Copyright (c) 2016 will freyman. All rights reserved.
//

#include "Mntr_StochasticBranchStateTimes.h"

#include <cstddef>
#include <cstdint>
#include <string>

#include "ArgumentRule.h"
#include "IntegerPos.h"
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
#include "StochasticBranchStateTimesMonitor.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore { class Tree; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }


using namespace RevLanguage;


Mntr_StochasticBranchStateTimes::Mntr_StochasticBranchStateTimes(void) : FileMonitor()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Mntr_StochasticBranchStateTimes* Mntr_StochasticBranchStateTimes::clone(void) const
{
    
    return new Mntr_StochasticBranchStateTimes(*this);
}


void Mntr_StochasticBranchStateTimes::constructInternalObject( void )
{
    
    const std::string& file_name      = static_cast<const RlString  &>( filename->getRevObject()           ).getValue();
    const std::string& sep            = static_cast<const RlString  &>( separator->getRevObject()          ).getValue();
    unsigned int       print_gen      = (int)static_cast<const IntegerPos   &>( printgen->getRevObject()      ).getValue();
    bool               app            = static_cast<const RlBoolean &>( append->getRevObject()             ).getValue();
    bool               wv             = static_cast<const RlBoolean &>( version->getRevObject()            ).getValue();
    
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* cdbdp_tdn = static_cast<const RevLanguage::Tree&>( cdbdp->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::Tree>* cdbdp_sn  = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree>* >( cdbdp_tdn );
    
    RevBayesCore::StateDependentSpeciationExtinctionProcess *sse_process = NULL;
    sse_process = dynamic_cast<RevBayesCore::StateDependentSpeciationExtinctionProcess*>( &cdbdp_sn->getDistribution() );
    sse_process->setSampleCharacterHistory( true );
    
    RevBayesCore::StochasticBranchStateTimesMonitor *m = new RevBayesCore::StochasticBranchStateTimesMonitor( cdbdp_sn, (std::uint64_t)print_gen, file_name, sep );
    m->setAppend( app );
    m->setPrintVersion( wv );
    
    delete value;
    value = m;
    
}


/** Get Rev type of object */
const std::string& Mntr_StochasticBranchStateTimes::getClassType(void)
{
    
    static std::string revType = "Mntr_StochasticBranchStateTimes";
    
    return revType;
    
}


/** Get class type spec describing type of object */
const TypeSpec& Mntr_StochasticBranchStateTimes::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Monitor::getClassTypeSpec() ) );
    
    return rev_type_spec;
    
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Mntr_StochasticBranchStateTimes::getMonitorName( void ) const
{
    
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "StochasticBranchStateTimes";
    
    return c_name;
    
}


/** Return member rules (no members) */
const MemberRules& Mntr_StochasticBranchStateTimes::getParameterRules(void) const
{
    
    static MemberRules monitor_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        monitor_rules.push_back( new ArgumentRule("cdbdp"          , TimeTree::getClassTypeSpec(),  "The character dependent birth-death process to monitor.",                      ArgumentRule::BY_REFERENCE, ArgumentRule::ANY ) );

        // add the rules from the base class
        const MemberRules &parentRules = FileMonitor::getParameterRules();
        monitor_rules.insert(monitor_rules.end(), parentRules.begin(), parentRules.end());

        rules_set = true;
    }
    
    return monitor_rules;
    
}


/** Get type spec */
const TypeSpec& Mntr_StochasticBranchStateTimes::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
    
}


/** Get type spec */
void Mntr_StochasticBranchStateTimes::printValue(std::ostream &o) const
{
    
    o << "Mntr_StochasticBranchStateTimes";
    
}


/** Set a member variable */
void Mntr_StochasticBranchStateTimes::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "cdbdp" )
    {
        cdbdp = var;
    }
    else
    {
        FileMonitor::setConstParameter(name, var);
    }
    
}
