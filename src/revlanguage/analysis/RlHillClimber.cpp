#include <cstddef>
#include <ostream>
#include <vector>

#include "ArgumentRules.h"
#include "HillClimber.h"
#include "RevObject.h"
#include "RlModel.h"
#include "RlMonitor.h"
#include "RlMove.h"
#include "RlHillClimber.h"
#include "RlString.h"
#include "TypeSpec.h"
#include "WorkspaceVector.h"
#include "MaximumLikelihoodAnalysis.h"
#include "Monitor.h"
#include "Move.h"
#include "RbHelpReference.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlMaximumLikelihoodAnalysis.h"

namespace RevBayesCore { class Model; }


using namespace RevLanguage;

HillClimber::HillClimber() : MaximumLikelihoodAnalysis()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
HillClimber* HillClimber::clone(void) const
{
    
    return new HillClimber(*this);
}


void HillClimber::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    const RevBayesCore::Model&                              mdl     = static_cast<const Model &>( model->getRevObject() ).getValue();
    const WorkspaceVector<Monitor> &                        ws_vec_mntr = static_cast<const WorkspaceVector<Monitor> &>( monitors->getRevObject() );
    RevBayesCore::RbVector<RevBayesCore::Monitor>           mntr;
    for ( size_t i = 0; i < ws_vec_mntr.size(); ++i )
    {
        mntr.push_back( ws_vec_mntr[i].getValue() );
    }
    const WorkspaceVector<Move>&                            ws_vec_mvs = static_cast<const WorkspaceVector<Move> &>( moves->getRevObject() );
    RevBayesCore::RbVector<RevBayesCore::Move>              mvs;
    for ( size_t i = 0; i < ws_vec_mvs.size(); ++i )
    {
        mvs.push_back( ws_vec_mvs[i].getValue() );
    }
    const std::string &                                     sched   = static_cast<const RlString &>( moveschedule->getRevObject() ).getValue();
    
    RevBayesCore::HillClimber *m = new RevBayesCore::HillClimber(mdl, mvs, mntr);
    m->setScheduleType( sched );
    
    value = new RevBayesCore::MaximumLikelihoodAnalysis(m);
    
}


/**
 * Get Rev type of object
 */
const std::string& HillClimber::getClassType(void)
{
    
    static std::string rev_type = "HillClimber";
    
    return rev_type;
}

/**
 * Get class type spec describing type of object
 */
const TypeSpec& HillClimber::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( MaximumLikelihoodAnalysis::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string HillClimber::getConstructorFunctionName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "HillClimber";
    
    return c_name;
}





/** Return member rules (no members) */
const MemberRules& HillClimber::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        // add the rules from the base class
        const MemberRules &parentRules = MaximumLikelihoodAnalysis::getParameterRules();
        memberRules.insert(memberRules.end(), parentRules.begin(), parentRules.end());
        
        rules_set = true;
    }
    
    return memberRules;
}


/** Get type spec */
const TypeSpec& HillClimber::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void HillClimber::printValue(std::ostream &o) const
{
    
    o << "Hill-Climber";
}


/** Set a member variable */
void HillClimber::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    MaximumLikelihoodAnalysis::setConstParameter(name, var);
    
}
