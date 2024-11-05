#include <cstddef>
#include <algorithm>
#include <functional>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ModelVector.h"
#include "Natural.h"
#include "Mcmcmc.h"
#include "OptionRule.h"
#include "Probability.h"
#include "RevObject.h"
#include "RbException.h"
#include "RealPos.h"
#include "RevNullObject.h"
#include "RlModel.h"
#include "RlMonitor.h"
#include "RlMove.h"
#include "RlMcmcmc.h"
#include "RlString.h"
#include "TypeSpec.h"
#include "WorkspaceVector.h"
#include "ModelObject.h"
#include "Monitor.h"
#include "MonteCarloAnalysis.h"
#include "MonteCarloAnalysisOptions.h"
#include "Move.h"
#include "RbBoolean.h"
#include "RbHelpReference.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlBoolean.h"
#include "RlMonteCarloAnalysis.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class Model; }


using namespace RevLanguage;

Mcmcmc::Mcmcmc() : MonteCarloAnalysis()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Mcmcmc* Mcmcmc::clone(void) const
{
    
    return new Mcmcmc(*this);
}


void Mcmcmc::constructInternalObject( void )
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
    int                                                     nchains = (int)static_cast<const Natural &>( num_chains->getRevObject() ).getValue();
    int                                                     si      = (int)static_cast<const Natural &>( swap_interval->getRevObject() ).getValue();
    double                                                  delta   = static_cast<const RealPos &>( delta_heat->getRevObject() ).getValue();
    int                                                     nreps   = (int)static_cast<const Natural &>( num_runs->getRevObject() ).getValue();
    const std::string &                                     comb    = static_cast<const RlString &>( combine_traces->getRevObject() ).getValue();
    int                                                     ntries  = (int)static_cast<const Natural &>( num_init_attempts->getRevObject() ).getValue();
    
    bool                                                    th      = static_cast<const RlBoolean &>( tune_heat->getRevObject() ).getValue();
    double                                                  tht     = static_cast<const Probability &>( tune_heat_target->getRevObject() ).getValue();
    const std::string &                                     sm      = static_cast<const RlString &>( swap_method->getRevObject() ).getValue();
    const std::string &                                     smo     = static_cast<const RlString &>( swap_mode->getRevObject() ).getValue();
    
    RevBayesCore::Mcmcmc *m = new RevBayesCore::Mcmcmc(mdl, mvs, mntr, sched, nchains, si, delta, ntries, th, tht, sm, smo);
    
    if (heat_temps->getRevObject() != RevNullObject::getInstance())
    {
        std::vector<double>                                 ht      = static_cast<const ModelVector<RealPos> &>( heat_temps->getRevObject() ).getValue();
        
        if(ht.size() != nchains)
        {
            throw RbException( "The number of heats provided must match the number of chains." );
        }
        
        std::sort(ht.begin(), ht.end(), std::greater<double>());
        
        if (ht[0] != 1.0)
        {
            throw RbException( "The highest heat will be used for the cold chain so it must be 1.0." );
        }
        
        m->setHeatsInitial(ht);
    }
    
    if (sm == "both")
    {
        std::int64_t si2 = 0;
        if (swap_interval2->getRevObject() != RevNullObject::getInstance())
        {
            si2 = static_cast<const Natural &>( swap_interval2->getRevObject() ).getValue();
        }
        else
        {
            si2 = si;
        }
        
        m->setSwapInterval2(si2);
    }
    
    
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
    
    value = new RevBayesCore::MonteCarloAnalysis(m,nreps,ct);
    
}


/**
 * Get Rev type of object
 */
const std::string& Mcmcmc::getClassType(void)
{
    
    static std::string rev_type = "Mcmcmc";
    
    return rev_type;
}

/**
 * Get class type spec describing type of object
 */
const TypeSpec& Mcmcmc::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( MonteCarloAnalysis::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Mcmcmc::getConstructorFunctionName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "mcmcmc";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Mcmcmc::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        // add the rules from the base class
        const MemberRules &parentRules = MonteCarloAnalysis::getParameterRules();
        memberRules.insert(memberRules.end(), parentRules.begin(), parentRules.end());
        
        memberRules.push_back( new ArgumentRule("nchains"    , Natural::getClassTypeSpec()  , "The number of chains to run.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural( long(4) ) ) );
        memberRules.push_back( new ArgumentRule("swapInterval" , Natural::getClassTypeSpec(), "The interval at which swaps (between neighbor chains if the swapMethod is neighbor or both, or between chains chosen randomly if the swapMethod is random) will be attempted.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural( long(10) )) );
        memberRules.push_back( new ArgumentRule("deltaHeat"    , RealPos::getClassTypeSpec(), "The delta parameter for the heat function.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RealPos(0.2) ) );
        memberRules.push_back( new ArgumentRule("heats"   , ModelVector<Probability>::getClassTypeSpec(), "The heats of chains, starting from the cold chain to hotter chains so the first value must be 1.0. If heats are specified directly then the delta parameter would be ignored.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, NULL ) );
        memberRules.push_back( new ArgumentRule("tuneHeat"  , RlBoolean::getClassTypeSpec() , "Should we tune the heats during burnin?", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new RlBoolean( false ) ) );
        memberRules.push_back( new ArgumentRule("tuneHeatTarget", Probability::getClassTypeSpec(), "The acceptance probability of adjacent chain swaps targeted by heats auto-tuning.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Probability( 0.23 ) ) );
        
        std::vector<std::string> options_swapMethod;
        options_swapMethod.push_back( "neighbor" );
        options_swapMethod.push_back( "random" );
        options_swapMethod.push_back( "both" );
        memberRules.push_back( new OptionRule( "swapMethod", new RlString("neighbor"), options_swapMethod, "The method used to swap chains." ) );
        
        memberRules.push_back( new ArgumentRule("swapInterval2" , Natural::getClassTypeSpec(), "The interval at which swaps between randomly chosen chains will be attempted (if the swapMethod is specified as both; otherwise it would be the same as swapInterval if not provided).", ArgumentRule::BY_VALUE, ArgumentRule::ANY, NULL) );
        
        std::vector<std::string> options_swapMode;
        options_swapMode.push_back( "single" );
        options_swapMode.push_back( "multiple" );
        memberRules.push_back( new OptionRule( "swapMode", new RlString("single"), options_swapMode, "Whether make a single attempt per swap interval or attempt multiple (= nchains-1 or choose(nchains,2) for neighbor or random swaps, respectively) times." ) );
        
        rules_set = true;
    }
    
    return memberRules;
}


/** Get type spec */
const TypeSpec& Mcmcmc::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void Mcmcmc::printValue(std::ostream &o) const
{
    
    o << "MetropolisCoupledMarkovChainMonteCarloAnalysis";
}


/** Set a member variable */
void Mcmcmc::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "nchains" )
    {
        num_chains = var;
    }
    else if ( name == "swapInterval" )
    {
        swap_interval = var;
    }
    else if ( name == "swapInterval2" )
    {
        swap_interval2 = var;
    }
    else if ( name == "deltaHeat" )
    {
        delta_heat = var;
    }
    else if ( name == "heats" )
    {
        heat_temps = var;
    }
    else if ( name == "tuneHeat" )
    {
        tune_heat = var;
    }
    else if ( name == "tuneHeatTarget" )
    {
        tune_heat_target = var;
    }
    else if ( name == "swapMethod" )
    {
        swap_method = var;
    }
    else if ( name == "swapMode" )
    {
        swap_mode = var;
    }
    else
    {
        MonteCarloAnalysis::setConstParameter(name, var);
    }
}
