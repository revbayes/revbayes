#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MonteCarloAnalysisOptions.h"
#include "Natural.h"
#include "Mcmc.h"
#include "RevObject.h"
#include "RealPos.h"
#include "RlModel.h"
#include "RlMonitor.h"
#include "RlMove.h"
#include "RlMcmc.h"
#include "RlString.h"
#include "TypeSpec.h"
#include "WorkspaceVector.h"
#include "Monitor.h"
#include "MonteCarloAnalysis.h"
#include "Move.h"
#include "RbHelpReference.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlMonteCarloAnalysis.h"

namespace RevBayesCore { class Model; }


using namespace RevLanguage;

Mcmc::Mcmc() : MonteCarloAnalysis()
{
        
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Mcmc* Mcmc::clone(void) const
{
    
    return new Mcmc(*this);
}


void Mcmc::constructInternalObject( void )
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
    int                                                     nreps   = (int)static_cast<const Natural &>( num_runs->getRevObject() ).getValue();
    int                                                     ntries  = (int)static_cast<const Natural &>( num_init_attempts->getRevObject() ).getValue();
    const std::string &                                     comb    = static_cast<const RlString &>( combine_traces->getRevObject() ).getValue();

    RevBayesCore::MonteCarloAnalysisOptions::TraceCombinationTypes ct = RevBayesCore::MonteCarloAnalysisOptions::SEQUENTIAL;
    if ( comb == "sequential" )
    {
        ct = RevBayesCore::MonteCarloAnalysisOptions::SEQUENTIAL;
    }
    else if ( comb == "mixed" )
    {
        ct = RevBayesCore::MonteCarloAnalysisOptions::MIXED;
    }
    else if ( comb == "none" )
    {
        ct = RevBayesCore::MonteCarloAnalysisOptions::NONE;
    }
    
    RevBayesCore::Mcmc *m = new RevBayesCore::Mcmc(mdl, mvs, mntr, ntries);
    m->setScheduleType( sched );
    
    double                                                  lHeat   = static_cast<const RealPos &>( likelihood_heat->getRevObject() ).getValue();
    double                                                  pHeat   = static_cast<const RealPos &>( posterior_heat->getRevObject() ).getValue();
    double                                                  prHeat  = static_cast<const RealPos &>( prior_heat->getRevObject() ).getValue();
    
    if (lHeat != 1.0)
    {
        m->setLikelihoodHeat(lHeat);
    }
    if (pHeat != 1.0)
    {
        m->setChainPosteriorHeat(pHeat);
    }
    if (prHeat != 1.0)
    {
        m->setChainPriorHeat(prHeat);
    }
    
    value = new RevBayesCore::MonteCarloAnalysis(m,nreps,ct);
    
}


/**
 * Get Rev type of object
 */
const std::string& Mcmc::getClassType(void)
{
    
    static std::string rev_type = "Mcmc";
    
    return rev_type;
}

/** 
 * Get class type spec describing type of object
 */
const TypeSpec& Mcmc::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( MonteCarloAnalysis::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Mcmc::getConstructorFunctionName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "mcmc";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Mcmc::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        // add the rules from the base class
        const MemberRules &parentRules = MonteCarloAnalysis::getParameterRules();
        memberRules.insert(memberRules.end(), parentRules.begin(), parentRules.end());
        
        memberRules.push_back( new ArgumentRule("priorHeat", RealPos::getClassTypeSpec(), "The power that the prior will be raised to.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RealPos(1.0) ) );
        memberRules.push_back( new ArgumentRule("likelihoodHeat", RealPos::getClassTypeSpec(), "The power that the likelihood will be raised to.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RealPos(1.0) ) );
        memberRules.push_back( new ArgumentRule("posteriorHeat", RealPos::getClassTypeSpec(), "The power that the posterior will be raised to.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RealPos(1.0) ) );
        
        rules_set = true;
    }
    
    return memberRules;
}


/** Get type spec */
const TypeSpec& Mcmc::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void Mcmc::printValue(std::ostream &o) const
{
    
    o << "MarkovChainMonteCarloAnalysis";
}


/** Set a member variable */
void Mcmc::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    if ( name == "likelihoodHeat" )
    {
        likelihood_heat = var;
    }
    else if ( name == "posteriorHeat" )
    {
        posterior_heat = var;
    }
    else if ( name == "priorHeat" )
    {
        prior_heat = var;
    }
    else
    {
        MonteCarloAnalysis::setConstParameter(name, var);
    }
    
}
