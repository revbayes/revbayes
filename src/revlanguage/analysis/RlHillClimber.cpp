#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ModelVector.h"
#include "Model.h"
#include "Natural.h"
#include "HillClimber.h"
#include "OptionRule.h"
#include "RevObject.h"
#include "RbException.h"
#include "Real.h"
#include "RealPos.h"
#include "RevNullObject.h"
#include "RlModel.h"
#include "RlMonitor.h"
#include "RlMove.h"
#include "RlHillClimber.h"
#include "RlString.h"
#include "TypeSpec.h"
#include "WorkspaceVector.h"


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


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> HillClimber::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Sebastian Hoehna" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string HillClimber::getHelpDescription(void) const
{
    std::string description = "";
    description += "The HillClimber analysis object keeps a model and the associated moves and monitors.";
    description += " The object is used to run Markov chain Monte Carlo (HillClimber) simulation on";
    description += " the model, using the provided moves, to obtain a sample of the posterior probability";
    description += " distribution. During the analysis, the monitors are responsible for sampling model parameters of interest.";
    
    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string HillClimber::getHelpDetails(void) const
{
    std::string details = "";
    details += " The HillClimber analysis object produced by a call to this function keeps copies of the model";
    details += " and the associated moves and monitors. The HillClimber analysis object is used to run Markov";
    details += " chain Monte Carlo (HillClimber) simulation on the model, using the provided moves,";
    details += " to obtain a sample of the posterior probability distribution. During the analysis,";
    details += " the monitors are responsible for sampling model parameters of interest.";
    
    return details;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string HillClimber::getHelpExample(void) const
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
    example += "# Create an HillClimber object\n";
    example += "myHillClimberObject = HillClimber( mymodel, monitors, moves)\n";
    example += "\n";
    example += "# Run a short analysis\n";
    example += "myHillClimberObject.burnin( generations = 400, tuningInterval = 100)\n";
    example += "myHillClimberObject.run( generations = 400)\n";
    example += "\n";
    example += "# print the summary of the operators (now tuned)\n";
    example += "myHillClimberObject.operatorSummary()\n";
    
    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> HillClimber::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
//    RevBayesCore::RbHelpReference ref1;
//    ref1.setCitation("Metropolis N, AW Rosenbluth, MN Rosenbluth, AH Teller, E Teller (1953). Equation of state calculations by fast computing machines. Journal of Chemical Physics, 21:1087-1092.");
//    references.push_back( ref1 );
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> HillClimber::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "SimulatedAnnealing" );
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string HillClimber::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Hill-Climber analysis object";
    
    return title;
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
