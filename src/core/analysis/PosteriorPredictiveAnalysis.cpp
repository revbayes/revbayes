#include <cstddef>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "DagNode.h"
#include "MaxIterationStoppingRule.h"
#include "MonteCarloAnalysis.h"
#include "MonteCarloSampler.h"
#include "MpiUtilities.h"
#include "PosteriorPredictiveAnalysis.h"
#include "RbException.h"
#include "Cloneable.h"
#include "Model.h"
#include "Parallelizable.h"
#include "RbFileManager.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "StoppingRule.h"

#ifdef RB_MPI
#include <mpi.h>
#endif


using namespace RevBayesCore;

PosteriorPredictiveAnalysis::PosteriorPredictiveAnalysis( const MonteCarloAnalysis &m, const path &fn ) : Cloneable( ), Parallelizable(),
    directory( fn ),
    processors_per_likelihood( 1 ),
    template_sampler( m )
{
    
//    // MPI settings
//    size_t active_proc = size_t( floor( pid / double(processors_per_likelihood)) );
//    m.setActivePID( active_proc );
//    m.setNumberOfProcesses( processors_per_likelihood );
    
    
}


PosteriorPredictiveAnalysis::PosteriorPredictiveAnalysis(const PosteriorPredictiveAnalysis &a) : Cloneable( a ), Parallelizable(a),
    directory( a.directory ),
    processors_per_likelihood( a.processors_per_likelihood ),
    template_sampler( a.template_sampler )
{
    
}


PosteriorPredictiveAnalysis::~PosteriorPredictiveAnalysis(void)
{
    
}



PosteriorPredictiveAnalysis* PosteriorPredictiveAnalysis::clone( void ) const
{
    
    return new PosteriorPredictiveAnalysis( *this );
}


void PosteriorPredictiveAnalysis::runAll(size_t gen)
{
    
    // print some information to the screen but only if we are the active process
    if ( process_active == true )
    {
        std::cout << std::endl;
        std::cout << "Running posterior predictive analysis ..." << std::endl;
    }
    
    // create the directory if necessary
    if ( not is_directory( directory) )
    {
        std::string errorStr = "";
        formatError(directory, errorStr);
        throw RbException()<<"Could not find file or path with name "<<directory<<"";
    }
    
    // set up a vector of strings containing the name or names of the files to be read
    std::vector<path> dir_names;
    if ( is_directory(directory) )
    {
        setStringWithNamesOfFilesInDirectory( directory, dir_names, false );
    }
    else
    {
        throw RbException()<<directory<<" is not a directory.";
    }

    size_t num_runs = dir_names.size();
    processors_per_likelihood = ceil( double(num_processes) / num_runs );
    size_t run_pid_start =  floor(  pid    / double(num_processes) * num_runs );
    size_t run_pid_end   =  floor( (pid+1) / double(num_processes) * num_runs );

    if ( run_pid_start == run_pid_end )
    {
        ++run_pid_end;
    }

    
//    int number_processes_per_run = int(run_pid_end) - int(run_pid_start) + 1;

    // set the processors for this analysis
    size_t active_proc = floor( pid / double(processors_per_likelihood) ) * processors_per_likelihood;
    template_sampler.setActivePID( active_proc, processors_per_likelihood );    
    
#ifdef RB_MPI
    MPI_Comm analysis_comm;
    MPI_Comm_split(MPI_COMM_WORLD, active_proc, pid, &analysis_comm);
#endif

    for ( size_t i = run_pid_start; i < run_pid_end; ++i)
    {
        
//        size_t run_pid_start = size_t(floor( double(i) / num_processes * num_runs ) );
//        size_t run_pid_end   = std::max( int(run_pid_start), int(floor( double(i+1) / num_processes * num_runs ) ) - 1);
        
        // create an independent copy of the analysis
        MonteCarloAnalysis *current_analysis = template_sampler.clone();

        // get the model of the analysis
        Model* current_model = current_analysis->getModel().clone();
        
        // get the DAG nodes of the model
        std::vector<DagNode*> &current_nodes = current_model->getDagNodes();
        
        // initialize values from files
        for (size_t j = 0; j < current_nodes.size(); ++j)
        {
            DagNode *the_node = current_nodes[j];
            if ( the_node->isClamped() == true )
            {
                the_node->setValueFromFile( dir_names[i] );
            }
                
        }

        // now set the model of the current analysis
#ifdef RB_MPI
        current_analysis->setModel( current_model, false, analysis_comm );
#else
        current_analysis->setModel( current_model, false );
#endif
        
        // disable the screen monitor
        current_analysis->disableScreenMonitors( true );

        // set the monitor index
        current_analysis->addFileMonitorExtension( dir_names[i].filename().string(), true);
                    
        // print some info
        if ( process_active == true )
        {
            size_t digits = size_t( ceil( log10( num_runs ) ) );
            std::cout << "Sim ";
            for (size_t d = size_t( ceil( log10( i+1.1 ) ) ); d < digits; ++d )
            {
                std::cout << " ";
            }
            std::cout << (i+1) << " / " << num_runs;
            std::cout << std::endl;
        }
        
        // run the i-th analysis
#ifdef RB_MPI
        runSim(current_analysis, gen, analysis_comm);
#else
        runSim(current_analysis, gen);
#endif
                
        // free memory
        delete current_analysis;
        
    }
    
    
#ifdef RB_MPI
    MPI_Comm_free(&analysis_comm);
    
    // to be safe, we should synchronize the random number generators
    MpiUtilities::synchronizeRNG( MPI_COMM_WORLD );
    
    // wait until all analysies are completed
    MPI_Barrier(MPI_COMM_WORLD);
#else
    MpiUtilities::synchronizeRNG(  );
#endif
    
}


#ifdef RB_MPI
void PosteriorPredictiveAnalysis::runSim(MonteCarloAnalysis *sampler, size_t gen, MPI_Comm &c)
#else
void PosteriorPredictiveAnalysis::runSim(MonteCarloAnalysis *sampler, size_t gen)
#endif
{
    
    // run the analysis
    RbVector<StoppingRule> rules;
    
    size_t currentGen = sampler->getCurrentGeneration();
    rules.push_back( MaxIterationStoppingRule(gen + currentGen) );
    
#ifdef RB_MPI
    sampler->run(gen, rules, c, 100, "", 0, false);
#else
    sampler->run(gen, rules, 100, "", 0, false);
#endif
    
}

