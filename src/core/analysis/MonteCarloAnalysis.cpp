#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <ostream>
#include <string>
#include <type_traits>
#include <vector>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "AbstractFileMonitor.h"
#include "DagNode.h"
#include "MonteCarloAnalysis.h"
#include "MonteCarloSampler.h"
#include "MpiUtilities.h"
#include "ProgressBar.h"
#include "RlUserInterface.h"
#include "Cloneable.h"
#include "Model.h"
#include "Monitor.h"
#include "MonteCarloAnalysisOptions.h"
#include "Parallelizable.h"
#include "RbException.h"
#include "RbFileManager.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "StoppingRule.h"
#include "Trace.h"


using namespace RevBayesCore;


/**
 * Constructor.
 *
 * \param[in]    m    The monte carlo sampler.
 */
MonteCarloAnalysis::MonteCarloAnalysis(MonteCarloSampler *m, size_t r, MonteCarloAnalysisOptions::TraceCombinationTypes tc) : Cloneable(), Parallelizable(),
    replicates( r ),
    runs(r,NULL),
    trace_combination( tc )
{
    
    runs[0] = m;

#ifdef RB_MPI
    MPI_Comm analysis_comm;
    MPI_Comm_split(MPI_COMM_WORLD, active_PID, pid, &analysis_comm);
    resetReplicates(analysis_comm);
#else
    resetReplicates();
#endif
    
}


MonteCarloAnalysis::MonteCarloAnalysis(const MonteCarloAnalysis &a) : Cloneable(), Parallelizable(a),
    replicates( a.replicates ),
    runs(a.replicates,NULL),
    trace_combination( a.trace_combination )
{
    
    // create replicate Monte Carlo samplers
    for (size_t i=0; i < replicates; ++i)
    {
        
        if ( a.runs[i] != NULL )
        {
            runs[i] = a.runs[i]->clone();
        }
        
    }
    
}




/**
 * Destructor. Nothing to do here
 */
MonteCarloAnalysis::~MonteCarloAnalysis(void)
{
    
    // free the runs
    for (size_t i = 0; i < replicates; ++i)
    {
        MonteCarloSampler *sampler = runs[i];
        delete sampler;
    }
    
}


/**
 * Overloaded assignment operator.
 * We need to keep track of the MonteCarloSamplers
 */
MonteCarloAnalysis& MonteCarloAnalysis::operator=(const MonteCarloAnalysis &a)
{
    Parallelizable::operator=(a);
    
    if ( this != &a )
    {
        
        // free the runs
        for (size_t i = 0; i < replicates; ++i)
        {
            MonteCarloSampler *sampler = runs[i];
            delete sampler;
        }
        runs = std::vector<MonteCarloSampler*>(a.replicates,NULL);
        
        replicates          = a.replicates;
        trace_combination   = a.trace_combination;
        
        // create replicate Monte Carlo samplers
        for (size_t i=0; i < replicates; ++i)
        {
            if ( a.runs[i] != NULL )
            {
                runs[i] = a.runs[i]->clone();
            }
            
        }
        
    }
    
    return *this;
}


/**
 * Set the model by delegating the model to the Monte Carlo samplers (replicates).
 */
void MonteCarloAnalysis::addFileMonitorExtension(const std::string &s, bool dir)
{
    
    // reset the counters for the move schedules
    for (size_t i=0; i<replicates; ++i)
    {
        runs[i]->addFileMonitorExtension(s, dir);
    }
    
}


/**
 * Add the monitors.
 */
void MonteCarloAnalysis::addMonitor(const Monitor &m)
{
    
    // remove the monitors for each replicate
    for (size_t i=0; i<replicates; ++i)
    {
        runs[i]->addMonitor( m );
    }
    
}

/** Run burnin and auto-tune */
#ifdef RB_MPI
void MonteCarloAnalysis::burnin(size_t generations, const MPI_Comm &analysis_comm, size_t tuningInterval, int verbose)
#else
void MonteCarloAnalysis::burnin(size_t generations, size_t tuningInterval, int verbose)
#endif
{
    
    // Initialize objects needed by chain
    for (size_t i=0; i<replicates; ++i)
    {
        
        if ( runs[i] != NULL )
        {
            runs[i]->initializeSampler();
        }
        
    }
    
    
    // reset the counters for the move schedules
    for (size_t i=0; i<replicates; ++i)
    {
        
        if ( runs[i] != NULL )
        {
            runs[i]->reset();
        }
        
    }
    
    // start the progress bar
    ProgressBar progress = ProgressBar(generations, 0);

    if ( verbose >= 1 && runs[0] != NULL && process_active == true )
    {
        // Let user know what we are doing
        std::stringstream ss;
        ss << "\n";
        ss << "Running burn-in phase of Monte Carlo sampler for " << generations << " iterations.\n";
        ss << "This simulation runs " << replicates << " independent replicate" << (replicates > 1 ? "s" : "") << ".\n";
        ss << runs[0]->getStrategyDescription();
        RBOUT( ss.str() );
        
        // Print progress bar (68 characters wide)
        progress.start();
    }
    
    
    // Run the chain
    for (size_t k=1; k<=generations; ++k)
    {
        
        if ( verbose >= 1 && process_active == true)
        {
            progress.update(k);
        }
        
        for (size_t i=0; i<replicates; ++i)
        {
            
            if ( runs[i] != NULL )
            {
                runs[i]->nextCycle(false);
                
                // check for autotuning
                if ( k % tuningInterval == 0 && k != generations )
                {
                    runs[i]->tune();
                }
                
            }
            
        }
        
    }
    
#ifdef RB_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    
    if ( verbose >= 1 && process_active == true )
    {
        progress.finish();
    }
    
#ifdef RB_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    
}



MonteCarloAnalysis* MonteCarloAnalysis::clone( void ) const
{
    
    return new MonteCarloAnalysis( *this );
}


void MonteCarloAnalysis::disableScreenMonitors(bool all)
{
    
    for (size_t i=0; i<replicates; ++i)
    {
        
        if ( runs[i] != NULL )
        {
            
            runs[i]->disableScreenMonitor(all, i);
        }
        
    }
    
}


size_t MonteCarloAnalysis::getCurrentGeneration( void ) const
{
    
    for (size_t i=0; i<replicates; ++i)
    {
        
        if ( runs[i] != NULL )
        {
            return runs[i]->getCurrentGeneration();
        }
        
    }
    
    return 0;
}


const Model& MonteCarloAnalysis::getModel( void ) const
{
    
    for (size_t i=0; i<replicates; ++i)
    {
        
        if ( runs[i] != NULL )
        {
            return runs[i]->getModel();
        }
        
    }
    
    return runs[0]->getModel();
}


void MonteCarloAnalysis::initializeFromCheckpoint(const path &checkpoint_file)
{
    
    for (size_t i = 0; i < replicates; ++i)
    {
        if ( runs[i] != NULL )
        {
            // first, set the checkpoint filename for the run
            if ( replicates > 1 && checkpoint_file != "" )
            {
                
                // create the run specific appendix
                std::stringstream ss;
                ss << "_run_" << (i+1);
                
                // assemble the new filename
                path run_checkpoint_file = appendToStem( checkpoint_file, ss.str() );
                
                // set the filename for the MCMC object
                runs[i]->setCheckpointFile( run_checkpoint_file );
            }
            else if ( not checkpoint_file.empty() )
            {
                // set the filename for the MCMC object
                runs[i]->setCheckpointFile( checkpoint_file );
            }
            
            // then, initialize the sample for that replicate
            runs[i]->initializeSamplerFromCheckpoint();
        }
    }
}


/**
 * Print out a summary of the current performance.
 */
void MonteCarloAnalysis::printPerformanceSummary( bool current_period ) const
{
    
#ifdef RB_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    
    if ( runs[0] != NULL )
    {
        runs[0]->printOperatorSummary( current_period );
    }
    
#ifdef RB_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    
}


/**
 * Remove the monitors.
 */
void MonteCarloAnalysis::removeMonitors( void )
{
    
    // remove the monitors for each replicate
    for (size_t i=0; i<replicates; ++i)
    {
        
        if ( runs[i] != NULL )
        {
            runs[i]->removeMonitors();
        }
        
    }
    
}


/**
 * Reset the replicates.
 */
#ifdef RB_MPI
void MonteCarloAnalysis::resetReplicates( const MPI_Comm &analysis_comm )
#else
void MonteCarloAnalysis::resetReplicates( void )
#endif
{

    // free the runs
    MonteCarloSampler *m = NULL;
    for (size_t i = 0; i < replicates; ++i)
    {
        MonteCarloSampler *sampler = runs[i];
        
        if ( m == NULL )
        {
            m = sampler;
        }
        
        if ( m != sampler )
        {
            delete sampler;
        }
        
        runs[i] = NULL;
        
    }
    
    
    if ( m == NULL )
    {
        throw RbException("Bug: No template sampler found!");
    }
    
    std::vector< size_t > replicate_indices_start = std::vector<size_t>(num_processes,0);
    std::vector< size_t > replicate_indices_end   = std::vector<size_t>(num_processes,0);
    
    for (size_t i=0; i<num_processes; ++i)
    {
        size_t this_replicate_start = size_t(floor( (double(i-active_PID) / num_processes ) * replicates ) );
        size_t this_replicate_end   = size_t(floor( (double(i+1-active_PID) / num_processes ) * replicates ) );

        replicate_indices_start[i] = this_replicate_start;
        replicate_indices_end[i]   = size_t( fmin( fmax(this_replicate_start+1,this_replicate_end), replicates ) );
        
    }
    
    // create replicate Monte Carlo samplers
    bool no_sampler_set = true;
    // making sure that initially only one core is used per sampler
    m->setActivePID( pid, 1 );
    for (size_t i = 0; i < replicates; ++i)
    {
        size_t replicate_pid_start = num_processes;
        size_t replicate_pid_end = 0;
        for (size_t j=0; j<num_processes; ++j)
        {

            if ( i >= replicate_indices_start[j] && i < replicate_indices_end[j] && replicate_pid_start > j )
            {
                replicate_pid_start = j;
            }
            if ( i >= replicate_indices_start[j] && i < replicate_indices_end[j] && replicate_pid_end < j )
            {
                replicate_pid_end = j;
            }
        }
        
        replicate_pid_start += active_PID;
        replicate_pid_end   += active_PID;

        int number_processes_per_replicate = int(replicate_pid_end) - int(replicate_pid_start) + 1;

        if ( pid >= replicate_pid_start && pid <= replicate_pid_end )
        {
            no_sampler_set = false;

            if ( i == 0 )
            {
                runs[i] = m;
//                runs[i]->setActivePID( replicate_pid_start, number_processes_per_replicate );
            }
            else
            {

                runs[i] = m->clone();
//                runs[i]->setActivePID( replicate_pid_start, number_processes_per_replicate );

            }
            
            runs[i]->setActivePID( replicate_pid_start, number_processes_per_replicate );
            //            runs[i]->setMasterSampler( i == 0 );
        }
        
    }
    if ( no_sampler_set == true )
    {
        runs[0] = m;
    }
    
    // disable the screen monitors for the replicates
    disableScreenMonitors( false );
    
    
    // we only need to tell the MonteCarloSamplers which replicate index they are if there is more than one replicate
    if ( replicates > 1 )
    {
        for (size_t i = 0; i < replicates; ++i)
        {
            
            if ( runs[i] != NULL )
            {
                
                std::stringstream ss;
                ss << "_run_" << (i+1);
                runs[i]->addFileMonitorExtension( ss.str(), false);
                
            }
            
        }
        
    }
    
    // get new random starting values
    size_t replicate_start = size_t(floor( (double(pid-active_PID) / num_processes ) * replicates ) ) + active_PID;
    
    RandomNumberGenerator *rng = GLOBAL_RNG;
    if(replicate_start > 0) rng->setSeed(rng->getSeed() + 2*replicate_start);    
    
    // redraw initial states for replicates
    for (size_t i = 0; i < replicates; ++i)
    {
        
        if ( i > 0 && runs[i] != NULL )
        {
            runs[i]->redrawStartingValues();
        }
        
        if ( runs[i] != NULL )
        {
            const std::vector<DagNode*> &this_nodes = runs[i]->getModel().getDagNodes();
            
            // touch all nodes
            for (size_t j=0; j<this_nodes.size(); ++j)
            {
                this_nodes[j]->touch();
            }
            
            // keep all nodes
            for (size_t j=0; j<this_nodes.size(); ++j)
            {
                this_nodes[j]->keep();
            }
        
        }
        
    }
    
    // to be safe, we should synchronize the random number generators
    // Sebastian: We cannot re-synchronize the RNG after we just shifted it.
    // If an analysis had all values preset, then each replicate would be identical!!!
//#ifdef RB_MPI
//    MpiUtilities::synchronizeRNG( analysis_comm );
//#else
//    MpiUtilities::synchronizeRNG(  );
//#endif
}



#ifdef RB_MPI
void MonteCarloAnalysis::run( size_t kIterations, RbVector<StoppingRule> rules, const MPI_Comm &analysis_comm, size_t tuning_interval, const path &checkpoint_file, size_t checkpoint_interval, int verbose )
#else
void MonteCarloAnalysis::run( size_t kIterations, RbVector<StoppingRule> rules, size_t tuning_interval, const path &checkpoint_file, size_t checkpoint_interval, int verbose )
#endif
{
    
    // get the current generation
    size_t gen = 0;
    for (size_t i=0; i<replicates; ++i)
    {
        
        if ( runs[i] != NULL )
        {
            gen = runs[i]->getCurrentGeneration();
            
            // also set the filename for checkpointing
            if ( replicates > 1 && checkpoint_file != "" )
            {
                
                // create the run specific appendix
                std::stringstream ss;
                ss << "_run_" << (i+1);
                
                // assemble the new filename
                auto run_checkpoint_file = appendToStem(checkpoint_file, ss.str());

                // set the filename for the MCMC object
                runs[i]->setCheckpointFile( run_checkpoint_file );
            }
            else if ( checkpoint_file != "" )
            {
                // set the filename for the MCMC object
                runs[i]->setCheckpointFile( checkpoint_file );
                
            }
            
        }
        
    }
    
    if ( process_active == true && runs[0] != NULL && verbose >= 1 )
    {
        // Let user know what we are doing
        std::stringstream ss;
        
        if ( runs[0]->getCurrentGeneration() == 0 )
        {
            ss << "\n";
            ss << "Running MCMC simulation\n";
        }
        else
        {
            ss << "Appending to previous MCMC simulation of " << runs[0]->getCurrentGeneration() << " iterations\n";
        }
        ss << "This simulation runs " << replicates << " independent replicate" << (replicates > 1 ? "s" : "") << ".\n";
        ss << runs[0]->getStrategyDescription();
        
        // Print the target values of stopping rules only if we have more than one or if the only one we have is not MaxIteration
        if (rules.size() > 1 or rules[0].printAsStatement(0, true) != "")
        {
            ss << "\n";
            ss << "Stopping rule" << (rules.size() > 1 ? "s" : "") << ":\n";
            for (size_t i=0; i<rules.size(); ++i)
            {
                std::string statement = rules[i].printAsStatement(0, true);
                if (statement == "")
                {
                    continue;
                }
                ss << "    " << statement;
            }
        }
        
        RBOUT( ss.str() );
    }
    
    // Start monitor(s)
    for (size_t i=0; i<replicates; ++i)
    {
        
        // Sebastian (2016/04/16): We should always reset the monitors so that the ETA starts fresh
        //        if ( runs[i] != NULL && runs[i]->getCurrentGeneration() == 0 )
        if ( runs[i] != NULL )
        {
            
            if ( i > 0 )
            {
                runs[i]->disableScreenMonitor(true, i);
            }
            
            runs[i]->startMonitors( kIterations, runs[i]->getCurrentGeneration() > 0 );
            
        }
        
    }
    
    // Sebastian: This is very important here!
    // We need to wait first for all processes and chains to have opened the filestreams
    // before we start printing (e.g., the headers) anything.
#ifdef RB_MPI
    // wait until all chains opened the monitor
    MPI_Barrier( analysis_comm );
#endif
    
    // reset the counters for the move schedules
    for (size_t i=0; i<replicates; ++i)
    {
        
        if ( runs[i] != NULL )
        {
            runs[i]->reset();
        }
        
    }
    
    // reset the stopping rules
    for (size_t i=0; i<rules.size(); ++i)
    {
        
        rules[i].setNumberOfRuns( replicates );
        rules[i].runStarted();
        
    }
    
    // If the conditions above are satisfied and we are resuming from a checkpoint, print the current stopping rule values too
    if ( runs[0]->getCurrentGeneration() != 0 and (rules.size() > 1 or rules[0].printAsStatement(0, true) != ""))
    {
        std::stringstream ss;
        ss << "Current value" << (rules.size() > 1 ? "s" : "") << ":\n";
        
        for (size_t i=0; i<rules.size(); ++i)
        {
            std::string to_parse = rules[i].printAsStatement( runs[0]->getCurrentGeneration() );
            // Delete the target info, which is redundant with respect to what we have already printed above
            std::string out = to_parse.substr(0, to_parse.find(" (target", 0));
            ss << "    " << out << "\n";
        }
        
        RBOUT( ss.str() );
    }
    
    // Write headers and print first line
    for (size_t i=0; i<replicates; ++i)
    {
        
        if ( runs[i] != NULL && runs[i]->getCurrentGeneration() == 0 )
        {
            
            runs[i]->writeMonitorHeaders( false );
            runs[i]->monitor(0);
            
        }
        else if ( runs[i] != NULL )
        {
            runs[i]->writeMonitorHeaders( runs[i]->getCurrentGeneration() > 0 );
        }
        
    }
    
    // Run the chain
    bool finished = false;
    bool converged = false;
    do {
        
        ++gen;
        for (size_t i=0; i<replicates; ++i)
        {            
            if ( runs[i] != NULL )
            {
                
                // @todo: #thread
                // This part should be done on several threads if possible
                // Sebastian: this call is very slow; a lot of work happens in nextCycle()
                runs[i]->nextCycle(true);
                
                // Monitor
                runs[i]->monitor(gen);
                
                // check for autotuning
                if ( tuning_interval != 0 && (gen % tuning_interval) == 0 )
                {                   
                    runs[i]->tune();                   
                }
                
                // check for checkpointing
                if ( checkpoint_interval != 0 && (gen % checkpoint_interval) == 0 )
                {                    
                    runs[i]->checkpoint();                    
                }             
            }           
        }
        
        converged = true;
        size_t numConvergenceRules = 0;
        
        // run the stopping test
        for (size_t i=0; i<rules.size(); ++i)
        {
            if ( rules[i].isConvergenceRule() )
            {
                converged &= rules[i].checkAtIteration(gen) && rules[i].stop(gen);
                ++numConvergenceRules;
            }
            else
            {
                if ( rules[i].checkAtIteration(gen) && rules[i].stop(gen) )
                {
                    finished = true;
                    break;
                }
            }
        }
        
        if (verbose > 1)
        {
            bool checkNow = false;
            
            for (size_t i=0; i<rules.size(); ++i)
            {
                if ( rules[i].isConvergenceRule() )
                {
                    // The non-convergence stopping rules (MaxTime and MaxIteration) are checked every single iteration.
                    // To avoid printing an enormous number of lines if these (either one of them or both) are the only rules we have,
                    // we will only print when at least one convergence rule wants us to.
                    checkNow |= rules[i].checkAtIteration(gen);
                }
            }
            
            if (checkNow)
            {
                std::stringstream ssConv;
                for (size_t i=0; i<rules.size(); ++i)
                {
                    // Prettify: insert a blank line before printing out the first stopping rule statement
                    if (i == 0)
                    {
                        ssConv << "\n";
                    }
                    ssConv << rules[i].printAsStatement(gen);
                }
                RBOUT( ssConv.str() );
            }
        }
        
        converged &= numConvergenceRules > 0;
        
    } while ( finished == false && converged == false);
    
#ifdef RB_MPI
    // wait until all replicates complete
    MPI_Barrier( analysis_comm );
#endif
    
    // Monitor
    for (size_t i=0; i<replicates; ++i)
    {
        
        if ( runs[i] != NULL )
        {
            runs[i]->finishMonitors( replicates, trace_combination );
        }
        
    }
    
    
#ifdef RB_MPI
    // wait until all replicates complete
    MPI_Barrier( analysis_comm );
    
    // to be safe, we should synchronize the random number generators
    MpiUtilities::synchronizeRNG( analysis_comm );
#else
    MpiUtilities::synchronizeRNG(  );
#endif
    
}


/**
 * Set the active PID of this specific Monte Carlo analysis.
 */
void MonteCarloAnalysis::setActivePIDSpecialized(size_t a, size_t n)
{
    
#ifdef RB_MPI
    MPI_Comm analysis_comm;
    MPI_Comm_split(MPI_COMM_WORLD, active_PID, pid, &analysis_comm);
    resetReplicates(analysis_comm);
#else
    resetReplicates();
#endif
    
}


/**
 * Set the model by delegating the model to the Monte Carlo samplers (replicates).
 */
#ifdef RB_MPI
void MonteCarloAnalysis::setModel(Model *m, bool redraw, const MPI_Comm &analysis_comm)
#else
void MonteCarloAnalysis::setModel(Model *m, bool redraw)
#endif
{
    
    // reset the counters for the move schedules
    for (size_t i=0; i<replicates; ++i)
    {
        if ( runs[i] != NULL )
        {
            
            if ( i == 0 )
            {
                runs[0]->setModel( m, redraw );
            }
            else
            {
                const Model *old_model = &runs[i]->getModel();
                delete old_model;
                
                runs[i]->setModel( m->clone(), redraw );
            }
            
        }
        
    }
    
    
#ifdef RB_MPI
    resetReplicates(analysis_comm);
#else
    resetReplicates();
#endif
    
}
