#include <cstddef>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "MonteCarloSampler.h"
#include "MoveSchedule.h"
#include "MpiUtilities.h"
#include "PowerPosteriorAnalysis.h"
#include "ProgressBar.h"
#include "RbException.h"
#include "RbFileManager.h"
#include "Cloneable.h"
#include "MonteCarloAnalysisOptions.h"
#include "Parallelizable.h"
#include "StringUtilities.h"


#ifdef RB_MPI
#include <mpi.h>
#endif

using namespace RevBayesCore;

PowerPosteriorAnalysis::PowerPosteriorAnalysis(MonteCarloSampler *m, const path &fn, size_t k) : Cloneable( ), Parallelizable(),
    filename( fn ),
    powers(),
    sampler( m ),
    sampleFreq( 100 ),
    processors_per_likelihood( k )
{
    
    initMPI();
}


PowerPosteriorAnalysis::PowerPosteriorAnalysis(const PowerPosteriorAnalysis &a) : Cloneable( a ), Parallelizable( a ),
    filename( a.filename ),
    powers( a.powers ),
    sampler( a.sampler->clone() ),
    sampleFreq( a.sampleFreq ),
    processors_per_likelihood( a.processors_per_likelihood )
{
    
}


PowerPosteriorAnalysis::~PowerPosteriorAnalysis(void)
{
    delete sampler;
}


/**
 * Overloaded assignment operator.
 * We need to keep track of the MonteCarloSamplers
 */
PowerPosteriorAnalysis& PowerPosteriorAnalysis::operator=(const PowerPosteriorAnalysis &a)
{
    Parallelizable::operator=( a );
    
    if ( this != &a )
    {
        
        // free the sampler
        delete sampler;
        
        filename                        = a.filename;
        powers                          = a.powers;
        sampler                         = a.sampler->clone();
        sampleFreq                      = a.sampleFreq;
        processors_per_likelihood       = a.processors_per_likelihood;
        
    }
    
    return *this;
}


/** Run burnin and autotune */
void PowerPosteriorAnalysis::burnin(size_t generations, size_t tuningInterval)
{
    
//    initMPI();
    
    // Initialize objects needed by chain
    sampler->initializeSampler();
    
    
    // reset the counters for the move schedules
    sampler->reset();
    
    // start the progress bar
    ProgressBar progress = ProgressBar(generations, 0);
    
    if ( process_active == true )
    {
        // Let user know what we are doing
        std::stringstream ss;
        ss << "\n";
        ss << "Running burn-in phase of Power Posterior sampler for " << generations << " iterations.\n";
        ss << sampler->getStrategyDescription();
        std::cout << ss.str() << std::endl;
    
        // Print progress bar (68 characters wide)
        progress.start();
    }
    
    
    // Run the chain
    for (size_t k=1; k<=generations; k++)
    {
        if ( process_active == true )
        {
            progress.update( k );
        }
        
        sampler->nextCycle(false);
            
        // check for autotuning
        if ( k % tuningInterval == 0 && k != generations )
        {
                
            sampler->tune();
            
        }
        
    }
    
    
    if ( process_active == true )
    {
        progress.finish();
    }
    
}



PowerPosteriorAnalysis* PowerPosteriorAnalysis::clone( void ) const
{
    
    return new PowerPosteriorAnalysis( *this );
}


void PowerPosteriorAnalysis::initMPI( void )
{
    
    size_t active_proc = size_t( floor( pid   / double(processors_per_likelihood)) ) * processors_per_likelihood;
    sampler->setActivePID( active_proc, processors_per_likelihood );
    
#ifdef RB_MPI
    // wait until all processes are complete
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    
}


void PowerPosteriorAnalysis::runAll(size_t gen, double burnin_fraction, size_t pre_burnin_generations, size_t tuning_interval)
{
    
    //    initMPI();
    
    if( gen < sampleFreq )
    {
        throw(RbException("Trying to run power posterior analysis for fewer generations than sampleFreq, no samples will be stored"));
    }
    
    // disable the screen monitor(s) if any
    sampler->disableScreenMonitor(true, 0);
    
    
    // print some information to the screen but only if we are the active process
    if ( process_active == true )
    {
        std::cout << std::endl;
        std::cout << "Running power posterior analysis ..." << std::endl;
    }
    
    // compute which block of the data this process needs to compute
    //    size_t stone_block_start = size_t(floor( (double(pid)   / num_processes ) * powers.size()) );
    //    size_t stone_block_end   = size_t(floor( (double(pid+1) / num_processes ) * powers.size()) );
    
    size_t stone_block_start =  floor( ( floor( pid   /double(processors_per_likelihood)) / (double(num_processes) / processors_per_likelihood) ) * powers.size() );
    size_t stone_block_end   =  floor( ( ceil( (pid+1)/double(processors_per_likelihood)) / (double(num_processes) / processors_per_likelihood) ) * powers.size() );
    
#ifdef RB_MPI
    // Wait until all processes are complete: this is to make sure we do not print step 1 before the message above
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    
    // Run the chain
    for (size_t i = stone_block_start; i < stone_block_end; ++i)
    {
    
        // run the i-th stone
        runStone(i, gen, burnin_fraction, pre_burnin_generations, tuning_interval);
        
    }
    
#ifdef RB_MPI
    // wait until all chains complete
    MPI_Barrier(MPI_COMM_WORLD);
    
    // to be safe, we should synchronize the random number generators
    MpiUtilities::synchronizeRNG( MPI_COMM_WORLD);
#else
    MpiUtilities::synchronizeRNG(  );
#endif
    
    if ( process_active == true )
    {
        summarizeStones();
    }
    
}



void PowerPosteriorAnalysis::runStone(size_t idx, size_t gen, double burnin_fraction, size_t pre_burnin_generations, size_t tuning_interval)
{
    // create the directory if necessary
    if (filename.filename().empty() or filename.filename_is_dot() or filename.filename_is_dot_dot())
    {
        throw RbException("Please provide a filename with an extension");
    }

    std::string stone_tag = "_stone_" + std::to_string(idx);

    path stoneFileName = appendToStem(filename, stone_tag);
    createDirectoryForFile( stoneFileName );
    
    std::ofstream outStream( stoneFileName.string() );
    outStream << "state\t" << "power\t" << "likelihood" << std::endl;

    // reset the sampler
    sampler->reset();

    size_t burnin = size_t( ceil( burnin_fraction*gen ) );
    size_t printInterval = size_t( round( fmax(1,gen/40.0) ) );
    
    /* Print output for users.
     * First, we will find the smallest PID such that the number of stones assigned to the corresponding process is equal to
     * ceil( powers.size() / num_processes ) rather than floor( powers.size() / num_processes ). This will be the process
     * that gets to print its status to the standard output. Note that process_active has a PID of 0, and does not always
     * satisfy this condition.
     */
    std::vector<size_t> start(num_processes);
    std::vector<size_t> end(num_processes);
    std::vector<size_t> ceil_pids;
    
    for (size_t i = 0; i < num_processes; ++i)
    {
        // see PowerPosteriorAnalysis::runAll()
        start[i] = floor( ( floor( i / double(processors_per_likelihood)) / (double(num_processes) / processors_per_likelihood) ) * powers.size() );
        end[i]   = floor( ( ceil((i+1)/double(processors_per_likelihood)) / (double(num_processes) / processors_per_likelihood) ) * powers.size() );
        if (end[i] - start[i] == ceil((double)powers.size() / num_processes))
        {
            ceil_pids.push_back(i);
        }
    }
    
    size_t pid_to_print = *std::min_element(ceil_pids.begin(), ceil_pids.end());
    
    if (pid == pid_to_print)
    {
        // Get the width of a given number in characters
        auto width = [&](size_t x) { return size_t( ceil( log10(x + 0.1) ) ); };
        
        if (num_processes == 1)
        {
            // Figure out how much whitespace the lines should be padded out with to keep everything neatly aligned
            size_t digits = width( powers.size() );
            
            std::cout << "Step ";
            for (size_t d = width(idx + 1); d < digits; d++)
            {
                std::cout << " ";
            }
            std::cout << (idx + 1);
        }
        else
        {
            // Get the values to print
            size_t step = idx + 1 - start[pid_to_print]; // make sure step counter starts from 1
            size_t lower_bound = (step - 1) * num_processes + 1;
            size_t upper_bound = std::min( step * num_processes, powers.size() );
            
            // Figure out how much whitespace the lines should be padded out with to keep everything neatly aligned
            size_t tmp = powers.size() - ceil((double)powers.size() / num_processes) + floor((double)powers.size() / num_processes);
            size_t max_upper_bound = (powers.size() % num_processes == 1) ? end[ end.size() - 2 ] : end[ end.size() - 1 ];
            size_t max_lower_bound = (powers.size() <= num_processes) ? 1 : ((powers.size() % num_processes > 1) ? tmp : start[ start.size() - 1 ]);
            size_t offset = width(max_lower_bound) + width(max_upper_bound) + 2; // add 2 for the "--"
            
            if (lower_bound != upper_bound)
            {
                std::cout << "Steps ";
                for (size_t d = width(lower_bound) + width(upper_bound) + 2; d < offset; d++)
                {
                    std::cout << " ";
                }
                std::cout << lower_bound << "--" << upper_bound;
            }
            else
            {
                std::cout << "Step ";
                for (size_t d = width(lower_bound); d < offset + 1; d++) // add 1 to account for "Step" vs "Steps"
                {
                    std::cout << " ";
                }
                std::cout << lower_bound;
            }
        }
    
        std::cout << " / " << powers.size() << "\t\t";
        std::cout.flush();
    }
    
    // set the power of this sampler
    sampler->setLikelihoodHeat( powers[idx] );
    
    sampler->addFileMonitorExtension( stone_tag, false);
    
    // let's do a pre-burnin
    for (size_t k=1; k<=pre_burnin_generations; k++)
    {
        
        sampler->nextCycle(false);
        
        // check for autotuning
        if ( k % tuning_interval == 0 && k != pre_burnin_generations )
        {
            sampler->tune();
        }
        
    }
    
    // Monitor
    sampler->startMonitors(gen, false);
    sampler->writeMonitorHeaders( false );
    sampler->monitor(0);
    
    double p = powers[idx];
    for (size_t k = 1; k <= gen; ++k)
    {
        
        if (pid == pid_to_print)
        {
            if (k % printInterval == 0)
            {
                std::cout << "*";
                std::cout.flush();
            }
        }
        
        sampler->nextCycle( true );

        // Monitor
        sampler->monitor(k);
        
        // sample the likelihood
        if ( k > burnin && k % sampleFreq == 0 )
        {
            // compute the joint likelihood
            double likelihood = sampler->getModelLnProbability(true);
            outStream << k << "\t" << p << "\t" << likelihood << std::endl;
        }
            
    }
    
    if (pid == pid_to_print)
    {
        std::cout << std::endl;
    }
    
    outStream.close();
    
    // Monitor
    sampler->finishMonitors( 1, MonteCarloAnalysisOptions::NONE );
    
}


void PowerPosteriorAnalysis::summarizeStones( void )
{
    // create the directory if necessary
    createDirectoryForFile( filename );
    
    std::ofstream outStream( filename.string() );
    outStream << "state\t" << "power\t" << "likelihood" << std::endl;

    // Append each stone
    for (size_t idx = 0; idx < powers.size(); ++idx)
    {
        std::string stone_tag = "_stone_" + std::to_string(idx);
        path stoneFileName = appendToStem( filename, stone_tag );

        // read the i-th stone
        std::ifstream inStream( stoneFileName.string() );
        if (inStream.is_open())
        {
            bool header = true;
            std::string line = "";
            while ( std::getline(inStream,line) )
            {
                // we need to skip the header line
                if ( header == true )
                {
                    header  = false;
                }
                else
                {
                    outStream << line << std::endl;
                }
            }
            inStream.close();
        }
        else
        {
            std::cerr << "Problem reading stone " << idx+1 << " from file " << stoneFileName << "." << std::endl;
        }

    }
    
    // closing the file stream
    outStream.close();

}


void PowerPosteriorAnalysis::setPowers(const std::vector<double> &p)
{
    powers = p;
}


void PowerPosteriorAnalysis::setSampleFreq(size_t sf)
{
    sampleFreq = sf;
}

