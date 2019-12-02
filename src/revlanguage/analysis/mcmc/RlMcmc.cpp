#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ModelVector.h"
#include "Model.h"
#include "MonteCarloAnalysisOptions.h"
#include "Natural.h"
#include "Mcmc.h"
#include "OptionRule.h"
#include "RevObject.h"
#include "RbException.h"
#include "Real.h"
#include "RealPos.h"
#include "RevNullObject.h"
#include "RlModel.h"
#include "RlMonitor.h"
#include "RlMove.h"
#include "RlMcmc.h"
#include "RlString.h"
#include "TypeSpec.h"
#include "WorkspaceVector.h"


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


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Mcmc::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Sebastian Hoehna" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Mcmc::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description;
    
    description += "The MCMC analysis object keeps a model and the associated moves and monitors.";
    description += " The object is used to run Markov chain Monte Carlo (MCMC) simulation on";
    description += " the model, using the provided moves, to obtain a sample of the posterior probability";
    description += " distribution. During the analysis, the monitors are responsible for sampling model parameters of interest.";
    
    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string Mcmc::getHelpDetails(void) const
{
    // create a variable for the description of the function
    std::string details = "";
    details += "The MCMC analysis object produced by a call to this function keeps copies of the model";
    details += " and the associated moves and monitors. The MCMC analysis object is used to run Markov";
    details += " chain Monte Carlo (MCMC) simulation on the model, using the provided moves,";
    details += " to obtain a sample of the posterior probability distribution. During the analysis,";
    details += " the monitors are responsible for sampling model parameters of interest.";

    return details;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Mcmc::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    
    example += "# Create a simple model (unclamped)\n";
    example += "a ~ exponential(1)\n";
    example += "mymodel = model(a)\n";
    example += "\n";
    example += "# Create a move vector and a monitor vector\n";
    example += "moves[1] = mvScale(a, lambda=1.0, weight=1.0)\n";
    example += "monitors[1] = mnFile(a,\"output/out.log\")\n";
    example += "\n";
    example += "# Create an mcmc object\n";
    example += "mymcmcObject = mcmc( mymodel, monitors, moves)\n";
    example += "\n";
    example += "# Run a short analysis\n";
    example += "mymcmcObject.burnin( generations = 400, tuningInterval = 100)\n";
    example += "mymcmcObject.run( generations = 400)\n";
    example += "\n";
    example += "# print the summary of the operators (now tuned)\n";
    example += "mymcmcObject.operatorSummary()\n";
    
    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> Mcmc::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    RevBayesCore::RbHelpReference ref1;
    ref1.setCitation("Metropolis N, AW Rosenbluth, MN Rosenbluth, AH Teller, E Teller (1953). Equation of state calculations by fast computing machines. Journal of Chemical Physics, 21:1087-1092.");
    references.push_back( ref1 );
    
    RevBayesCore::RbHelpReference ref2;
    ref2.setCitation("Hastings WK (1970) Monte Carlo sampling methods using Markov chains and their applications. Biometrika, 57:97-109.");
    references.push_back( ref2 );
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Mcmc::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "mcmcmc" );
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Mcmc::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "MCMC analysis object";
    
    return title;
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
