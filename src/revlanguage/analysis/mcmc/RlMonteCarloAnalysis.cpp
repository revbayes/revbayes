
#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MaxIterationStoppingRule.h"
#include "MonteCarloAnalysis.h"
#include "Natural.h"
#include "OptionRule.h"
#include "RbException.h"
#include "RlMonteCarloAnalysis.h"
#include "RlModel.h"
#include "RlModelTrace.h"
#include "RlMonitor.h"
#include "RlMove.h"
#include "RlStoppingRule.h"
#include "RlString.h"
#include "StoppingRule.h"
#include "TypeSpec.h"
#include "WorkspaceVector.h"
#include "Argument.h"
#include "MemberProcedure.h"
#include "MethodTable.h"
#include "RbBoolean.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "RevNullObject.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlBoolean.h"
#include "RlUtils.h"
#include "Trace.h"
#include "WorkspaceToCoreWrapperObject.h"
#include "RlUserInterface.h"

using namespace RevLanguage;

MonteCarloAnalysis::MonteCarloAnalysis(void) : WorkspaceToCoreWrapperObject<RevBayesCore::MonteCarloAnalysis>( )
{
    
    initializeMethods();
   
}

MonteCarloAnalysis::MonteCarloAnalysis(RevBayesCore::MonteCarloAnalysis *m) : WorkspaceToCoreWrapperObject<RevBayesCore::MonteCarloAnalysis>( m )
{
    
    initializeMethods();
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
MonteCarloAnalysis* MonteCarloAnalysis::clone(void) const
{
    
    return new MonteCarloAnalysis(*this);
}


//void MonteCarloAnalysis::constructInternalObject( void )
//{
//    // we free the memory first
//    delete value;
//    
//    // now allocate a new MonteCarloAnalysis object
//    const RevBayesCore::Model&                              mdl     = static_cast<const Model &>( model->getRevObject() ).getValue();
//    const RevBayesCore::RbVector<RevBayesCore::Move>&       mvs     = static_cast<const WorkspaceVector<Move> &>( moves->getRevObject() ).getVectorRbPointer();
//    const RevBayesCore::RbVector<RevBayesCore::Monitor>&    mntr    = static_cast<const WorkspaceVector<Monitor> &>( monitors->getRevObject() ).getVectorRbPointer();
//    const std::string &                                     sched   = static_cast<const RlString &>( moveSchedule->getRevObject() ).getValue();
//    value = new RevBayesCore::MonteCarloAnalysis(mdl, mvs, mntr);
//    value->setScheduleType( sched );
//}

void deprecateUnderPrior(const std::string& method, const Argument& arg)
{
    static int warned = 0;

    auto& underPriorObject = arg.getVariable()->getRevObject();

    // If the argument is NULL the the user didn't specify `underPrior=value`.
    if (underPriorObject == RevNullObject::getInstance()) return;

    // What value did the user specify?
    auto underPrior = static_cast<const RlBoolean &>( underPriorObject ).getValue();

    // Skip if this is a warning and we've already showed the warning twice.
    if (not underPrior and warned > 2) return;

    RBOUT("\nWarning! The `underPrior` argument to `"+method+"` was removed in version 1.2.5 and now has no effect.");
    RBOUT("         You can accomplish the same effect as `underPrior=TRUE` with `model.ignoreAllData()`.");
    RBOUT("         See help on `model` for more information.\n");

    if (underPrior)
        throw RbException()<<method<<": argument `underPrior=TRUE` has no effect.";

    warned++;
}

/* Map calls to member methods */
RevPtr<RevVariable> MonteCarloAnalysis::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    
    if (name == "run")
    {
        found = true;
        
        int currentGen = int(value->getCurrentGeneration());
        
        size_t args_index = 0;
        
        RevBayesCore::RbVector<RevBayesCore::StoppingRule> rules;
        int gen = (int)static_cast<const Natural &>( args[args_index++].getVariable()->getRevObject() ).getValue();
        if(gen > 0) {
            gen += currentGen;
            rules.push_back( RevBayesCore::MaxIterationStoppingRule(gen) );
        }
        
        // get the member with given index
        if ( args[args_index++].getVariable()->getRevObject() != RevNullObject::getInstance() )
        {
            const WorkspaceVector<StoppingRule>& ws_vec = static_cast<const WorkspaceVector<StoppingRule> &>( args[1].getVariable()->getRevObject() );
            for ( size_t i = 0; i < ws_vec.size(); ++i )
            {
                rules.push_back( ws_vec[i].getValue() );
            }
        }

        if(rules.size() < 1) {
            throw RbException("Analysis requires a maximum number of generations or at least one stopping rule");
        }
        
        // the tuning interval (0 by default)
        int tuning_interval = (int)static_cast<const Natural &>( args[args_index++].getVariable()->getRevObject() ).getValue();

        // the checkoint file (empty by default) and the checkoint interval (0 by default)
        const std::string checkpoint_file = static_cast<const RlString &>( args[args_index++].getVariable()->getRevObject() ).getValue();
        long checkpoint_interval = static_cast<const Natural &>( args[args_index++].getVariable()->getRevObject() ).getValue();
        if ( (checkpoint_file != "" && checkpoint_interval == 0) || (checkpoint_file == "" && checkpoint_interval != 0) )
        {
            throw RbException("For checkpointing you have to provide both the checkpoint file and the checkpoint frequency (interval).");
        }

        deprecateUnderPrior("run", args[args_index++]);

#ifdef RB_MPI
        value->run( gen, rules, MPI_COMM_WORLD, tuning_interval, checkpoint_file, checkpoint_interval );
#else
        value->run( gen, rules, tuning_interval, checkpoint_file, checkpoint_interval );
#endif
        return NULL;
    }
    else if (name == "burnin")
    {
        found = true;
        
        // get the member with give index
        int gen = (int)static_cast<const Natural &>( args[0].getVariable()->getRevObject() ).getValue();
        int tuningInterval = (int)static_cast<const Natural &>( args[1].getVariable()->getRevObject() ).getValue();

        deprecateUnderPrior("burnin", args[2]);

#ifdef RB_MPI
        value->burnin( gen, MPI_COMM_WORLD, tuningInterval );
#else
        value->burnin( gen, tuningInterval );
#endif
        
        return NULL;
    }
    else if ( name == "operatorSummary")
    {
        found = true;
        
        bool current_period = static_cast<const RlBoolean &>( args[0].getVariable()->getRevObject() ).getValue();
        
        value->printPerformanceSummary( current_period );
        
        return NULL;
    }
    else if ( name == "initializeFromCheckpoint")
    {
        found = true;
        
        const std::string &checkpoint_filename = static_cast<const RlString &>( args[0].getVariable()->getRevObject() ).getValue();
        
        value->initializeFromCheckpoint( checkpoint_filename );
        
        return NULL;
    }
    
    return RevObject::executeMethod( name, args, found );
}


/** Get Rev type of object */
const std::string& MonteCarloAnalysis::getClassType(void)
{
    
    static std::string rev_type = "MonteCarloAnalysis";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& MonteCarloAnalysis::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( WorkspaceToCoreWrapperObject<RevBayesCore::MonteCarloAnalysis>::getClassTypeSpec() ) );
    
    return rev_type_spec;
}



/** Return member rules (no members) */
const MemberRules& MonteCarloAnalysis::getParameterRules(void) const
{
    
    static MemberRules member_rules;
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        
        member_rules.push_back( new ArgumentRule("model"   , Model::getClassTypeSpec()                   , "The model graph.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        member_rules.push_back( new ArgumentRule("monitors", WorkspaceVector<Monitor>::getClassTypeSpec(), "The monitors used for this analysis.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        member_rules.push_back( new ArgumentRule("moves"   , WorkspaceVector<Move>::getClassTypeSpec()   , "The moves used for this analysis.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        
        std::vector<std::string> options_schedule;
        options_schedule.push_back( "sequential" );
        options_schedule.push_back( "random" );
        options_schedule.push_back( "single" );
        
        member_rules.push_back( new OptionRule( "moveschedule", new RlString( "random" ), options_schedule, "The strategy how the moves are used." ) );
        
        // the number of replicate runs (in parallel)
        member_rules.push_back( new ArgumentRule("nruns"   , Natural::getClassTypeSpec(), "The number of replicate analyses.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(1) ) );
        
        // the way how we combine several traces at the end of the run
        std::vector<std::string> options_combine;
        options_combine.push_back( "sequential" );
        options_combine.push_back( "mixed" );
        options_combine.push_back( "none" );
        member_rules.push_back( new OptionRule( "combine", new RlString( "none" ), options_combine, "How should we combine the traces once the simulation is finished." ) );

        // the number of tries to initialize the MCMC until it fails
        member_rules.push_back( new ArgumentRule("ntries"   , Natural::getClassTypeSpec(), "The number of initialization attempts.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(1000) ) );

        rules_set = true;
    }
    
    return member_rules;
}

/** Get type spec */
const TypeSpec& MonteCarloAnalysis::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


void MonteCarloAnalysis::initializeMethods()
{
    
    ArgumentRules* run_arg_rules = new ArgumentRules();
    run_arg_rules->push_back( new ArgumentRule( "generations", Natural::getClassTypeSpec(), "The number of generations to run.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(0L) ) );
    run_arg_rules->push_back( new ArgumentRule( "rules", WorkspaceVector<StoppingRule>::getClassTypeSpec(), "The rules when to automatically stop the run.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, NULL ) );
    run_arg_rules->push_back( new ArgumentRule( "tuningInterval", Natural::getClassTypeSpec(), "The interval when to update the tuning parameters of the moves.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(0L)  ) );
    run_arg_rules->push_back( new ArgumentRule( "checkpointFile", RlString::getClassTypeSpec(), "The filename for the checkpoint file.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString("")  ) );
    run_arg_rules->push_back( new ArgumentRule( "checkpointInterval", Natural::getClassTypeSpec(), "The interval when to write parameters values to a files for checkpointing.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(0L)  ) );
    run_arg_rules->push_back( new ArgumentRule( "underPrior" , RlBoolean::getClassTypeSpec(), "DEPRECATED", ArgumentRule::BY_VALUE, ArgumentRule::ANY, nullptr ) );
    methods.addFunction( new MemberProcedure( "run", RlUtils::Void, run_arg_rules) );

    ArgumentRules* burninArgRules = new ArgumentRules();
    burninArgRules->push_back( new ArgumentRule( "generations"   , Natural::getClassTypeSpec(), "The number of generation to run this burnin simulation.", ArgumentRule::BY_VALUE, ArgumentRule::ANY  ) );
    burninArgRules->push_back( new ArgumentRule( "tuningInterval", Natural::getClassTypeSpec(), "The interval when to update the tuning parameters of the moves.", ArgumentRule::BY_VALUE, ArgumentRule::ANY  ) );
    burninArgRules->push_back( new ArgumentRule( "underPrior" , RlBoolean::getClassTypeSpec(), "DEPRECATED", ArgumentRule::BY_VALUE, ArgumentRule::ANY, nullptr ) );
    methods.addFunction( new MemberProcedure( "burnin", RlUtils::Void, burninArgRules) );

    ArgumentRules* operatorSummaryArgRules = new ArgumentRules();
    operatorSummaryArgRules->push_back( new ArgumentRule( "currentPeriod" , RlBoolean::getClassTypeSpec(), "Should the operator summary (number of tries and acceptance, and the acceptance ratio) of only the current period (i.e., after the last tuning) be printed?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );
    methods.addFunction( new MemberProcedure( "operatorSummary", RlUtils::Void, operatorSummaryArgRules) );
    
    ArgumentRules* initialize_checkpoint_arg_rules = new ArgumentRules();
    initialize_checkpoint_arg_rules->push_back( new ArgumentRule("checkpointFile", RlString::getClassTypeSpec(), "The checkpoint filename.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberProcedure( "initializeFromCheckpoint", RlUtils::Void, initialize_checkpoint_arg_rules) );
    
}


/** 
 * Print value
 */
void MonteCarloAnalysis::printValue(std::ostream &o) const
{
    
    o << "MonteCarloAnalysis";
}


/** 
 * Set a member variable 
 */
void MonteCarloAnalysis::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "model")
    {
        model = var;
    }
    else if ( name == "moves")
    {
        moves = var;
    }
    else if ( name == "monitors")
    {
        monitors = var;
    }
    else if ( name == "moveschedule")
    {
        moveschedule = var;
    }
    else if ( name == "nruns")
    {
        num_runs = var;
    }
    else if ( name == "combine")
    {
        combine_traces = var;
    }
    else if ( name == "ntries")
    {
        num_init_attempts = var;
    }
    else
    {
        WorkspaceToCoreWrapperObject<RevBayesCore::MonteCarloAnalysis>::setConstParameter(name, var);
    }
    
}
