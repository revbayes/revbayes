#include <cstddef>
#include <cmath>
#include <iostream>
#include <vector>

#include "DagNode.h"
#include "MaximumLikelihoodAnalysis.h"
#include "BootstrapAnalysis.h"
#include "Cloneable.h"
#include "Model.h"
#include "Parallelizable.h"


using namespace RevBayesCore;

BootstrapAnalysis::BootstrapAnalysis( const MaximumLikelihoodAnalysis &m, size_t n ) : Cloneable( ), Parallelizable( ),
    num_runs( n )
{
    
    // remove all monitors if there are any
    mle_analysis = m.clone();
    
}


BootstrapAnalysis::BootstrapAnalysis(const BootstrapAnalysis &a) : Cloneable( a ), Parallelizable( a ),
    num_runs( a.num_runs )
{
    
    mle_analysis = a.mle_analysis->clone();
}


BootstrapAnalysis::~BootstrapAnalysis(void)
{
    // free the runs
    delete mle_analysis;
    
}


/**
 * Overloaded assignment operator.
 * We need to keep track of the MonteCarloSamplers
 */
BootstrapAnalysis& BootstrapAnalysis::operator=(const BootstrapAnalysis &a)
{
    Parallelizable::operator=( a );
    
    if ( this != &a )
    {
        
        // free the runs
        delete mle_analysis;

        num_runs                    = a.num_runs;
        
        
        // create replicate Monte Carlo samplers
        mle_analysis = a.mle_analysis->clone();
        
        
    }
    
    return *this;
}


BootstrapAnalysis* BootstrapAnalysis::clone( void ) const
{
    
    return new BootstrapAnalysis( *this );
}


void BootstrapAnalysis::runAll(double epsilon)
{
    // disable all the screen monitors
    mle_analysis->disableScreenMonitors( true );
    
    // start the monitors
    mle_analysis->startMonitors();
    
    // print some information to the screen but only if we are the active process
    if ( process_active )
    {
        std::cout << std::endl;
        std::cout << "Running Bootstrap analysis ..." << std::endl;
    }
    
    // compute which block of the runs this process needs to compute
    size_t run_block_start = size_t(floor( (double(pid)   / num_processes ) * num_runs) );
    size_t run_block_end   = size_t(floor( (double(pid+1) / num_processes ) * num_runs) );
    //    size_t stone_block_size  = stone_block_end - stone_block_start;
    
    // Run the chain
    for (size_t i = run_block_start; i < run_block_end; ++i)
    {
        
        // run the i-th stone
        runSim(i, epsilon);
        
    }
    
    
    // close the monitors
    mle_analysis->finishMonitors();
    
}



void BootstrapAnalysis::runSim(size_t idx, double epsilon)
{
    // print some info
    if ( process_active )
    {
        size_t digits = size_t( ceil( log10( num_runs ) ) );
        std::cout << "Sim ";
        for (size_t d = size_t( ceil( log10( idx+1.1 ) ) ); d < digits; d++ )
        {
            std::cout << " ";
        }
        std::cout << (idx+1) << " / " << num_runs;
        std::cout << "\t\t";
        
        std::cout << std::endl;
    }
    
    // get the model of the analysis
    Model& current_model = mle_analysis->getModel();
    
    // get the DAG nodes of the model
    std::vector<DagNode *> current_ordered_nodes = current_model.getOrderedStochasticNodes();
    
    for (size_t j = 0; j < current_ordered_nodes.size(); ++j)
    {
        DagNode *the_node = current_ordered_nodes[j];
        
        if ( the_node->isStochastic() == true && the_node->isClamped() )
        {
            the_node->bootstrap();
            the_node->touch();
        }
        
    }
    
    
    // now run the analysis
    mle_analysis->run(epsilon, false);
    
    // save the estimate
    mle_analysis->monitor( idx );
    
    
//    for (size_t j = 0; j < current_ordered_nodes.size(); ++j)
//    {
//        DagNode *the_node = current_ordered_nodes[j];
//        
//        if ( the_node->isStochastic() == true && the_node->isClamped() )
//        {
//            the_node->bootstrap();
//        }
//        
//    }
    
}

