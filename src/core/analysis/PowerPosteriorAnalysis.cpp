#include <algorithm>
#include <cstddef>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
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
#include "RlUserInterface.h"
#include "StringUtilities.h"

#ifdef _WIN32
#include <windows.h>
#endif

#ifdef RB_MPI
#include <mpi.h>
#endif

using namespace RevBayesCore;

PowerPosteriorAnalysis::PowerPosteriorAnalysis(MonteCarloSampler *m, const path &fn, size_t k) : Cloneable( ), Parallelizable(),
    filename( fn ),
    powers(),
    sampler( m ),
    sampleFreq( 100 ),
    processors_per_likelihood( k ),
    resume_from_checkpoint( false ),
    ckp_stone_file()
{
    
    initMPI();
}


PowerPosteriorAnalysis::PowerPosteriorAnalysis(const PowerPosteriorAnalysis &a) : Cloneable( a ), Parallelizable( a ),
    filename( a.filename ),
    powers( a.powers ),
    sampler( a.sampler->clone() ),
    sampleFreq( a.sampleFreq ),
    processors_per_likelihood( a.processors_per_likelihood ),
    resume_from_checkpoint( a.resume_from_checkpoint ),
    ckp_stone_file( a.ckp_stone_file )
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
        resume_from_checkpoint          = a.resume_from_checkpoint;
        ckp_stone_file                  = a.ckp_stone_file;
        
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
        if ( tuningInterval != 0 && (k % tuningInterval) == 0 && k != generations )
        {
            sampler->tune();
        }
    }
    
    
    if ( process_active == true )
    {
        progress.finish();
    }
    
}


void PowerPosteriorAnalysis::checkpoint( size_t stone_idx, const path &base_checkpoint_file_name, size_t planned_burnin ) const
{
    path stone_file_name = appendToStem( base_checkpoint_file_name, "_stone_" + std::to_string(stone_idx + 1) );
    
    sampler->setCheckpointFile( stone_file_name );
    sampler->checkpoint();
    
    // save the power and planned burnin for this stone
    path stone_info_file_name = appendToStem( stone_file_name, "_stone_info" );

    path tmp_stone_info_file_name = stone_info_file_name.parent_path() / ("." + stone_info_file_name.filename().string() + ".tmp");
    std::ofstream out_stream_mcmc( tmp_stone_info_file_name.string() );
    out_stream_mcmc << "power = " << powers[stone_idx] << std::endl;
    out_stream_mcmc << "planned_burnin = " << planned_burnin << std::endl;

    out_stream_mcmc.close();
    const bool ok = out_stream_mcmc.good();
    if ( !ok )
    {
        RBOUT( "Warning: failed to write \"" + stone_info_file_name.string() + "\"; keeping existing file." );
        std::error_code ec;
        std::filesystem::remove( tmp_stone_info_file_name, ec );
    }
    else
#ifdef _WIN32
    if ( MoveFileExW( tmp_stone_info_file_name.wstring().c_str(), stone_info_file_name.wstring().c_str(), MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH ) == 0 )
    {
        throw RbException() << "Could not replace checkpoint file " << stone_info_file_name;
    }
#else
    std::filesystem::rename(tmp_stone_info_file_name, stone_info_file_name);
#endif
}


PowerPosteriorAnalysis* PowerPosteriorAnalysis::clone( void ) const
{
    
    return new PowerPosteriorAnalysis( *this );
}


std::vector<double> PowerPosteriorAnalysis::getPowers( void ) const
{
    return powers;
}


/**
 * Initialize the analysis from checkpoint files. Note that this function is much "lighter" than the equivalent functions of the
 * Mcmc and Mcmcmc classes, as it does not itself load the sampler state: we need to leave that up to runStone() / runAll().
 * This is because PowerPosteriorAnalysis only has one MonteCarloSampler* for all stones, unlike the previous two classes that
 * have a different one for each run / chain. We therefore have no place to hold the full state of each stone. Iterating over
 * stone_indices would just load the sampler state for each stone, overwriting the previous one. Instead, this function just
 * (1) checks that we can find the checkpoint files for the requested stones, (2) records their full paths for runStone() / runAll()
 * to use, and (3) sets the resume_from_checkpoint variable that tells those two functions to not start fresh.
 *
 * @param base_checkpoint_file_name The base name of the checkpoint files.
 * @param stone_indices The 1-based indices of the stones to restore from checkpoint.
 */
void PowerPosteriorAnalysis::initializeFromCheckpoint(const path &base_checkpoint_file_name, const std::vector<size_t> &stone_indices )
{
    if ( stone_indices.empty() )
    {
        throw RbException() << "Specify which stones should be restored from checkpoint.";
    }
    
    // sort the indices and check for duplicates
    std::vector<size_t> sorted_indices( stone_indices.size() );
    std::sort( sorted_indices.begin(), sorted_indices.end() );
    sorted_indices.erase( std::unique( sorted_indices.begin(), sorted_indices.end() ), sorted_indices.end() );
    
    ckp_stone_file.clear();
    
    for (size_t idx : sorted_indices)
    {
        if ( idx >= powers.size() )
        {
            throw RbException() << "Stone index " << (idx + 1) << " is larger than the number of stones (" << powers.size() << ").";
        }
        
        path stone_file_name = appendToStem( base_checkpoint_file_name, "_stone_" + std::to_string(idx + 1) );
        if ( !is_regular_file( stone_file_name ) )
        {
            throw RbException() << "No checkpoint file found for stone " << (idx + 1) << ".";
        }
        
        ckp_stone_file[idx] = stone_file_name;
    }
    
    resume_from_checkpoint = true;
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


void PowerPosteriorAnalysis::runAll(size_t gen, double burnin_fraction, size_t pre_burnin_generations, size_t tuning_interval, const path &checkpoint_file, size_t checkpoint_interval)
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
        std::cout << "Stone pre-burnin: -, burnin: +, sampling phase: *" << std::endl;
    }
    
#ifdef RB_MPI
    // Wait until all processes are complete: this is to make sure we do not print step 1 before the message above
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    
    if ( resume_from_checkpoint )
    {
        std::vector<size_t> resurrection_indices;
        for (auto& [k, v] : ckp_stone_file) resurrection_indices.push_back(k);
        
        size_t m = resurrection_indices.size();

        size_t block_start = size_t( floor( ( floor( pid   / double(processors_per_likelihood)) / (double(num_processes) / processors_per_likelihood) ) * m ) );
        size_t block_end   = size_t( floor( ( ceil( (pid + 1) / double(processors_per_likelihood)) / (double(num_processes) / processors_per_likelihood) ) * m ) );
        
        for (size_t j = block_start; j < block_end; ++j)
        {
            runStone( resurrection_indices[j], gen, burnin_fraction, pre_burnin_generations, tuning_interval, false, checkpoint_file, checkpoint_interval );
        }
    }
    else
    {
        // compute which block of the data this process needs to compute
        //    size_t stone_block_start = size_t(floor( (double(pid)   / num_processes ) * powers.size()) );
        //    size_t stone_block_end   = size_t(floor( (double(pid+1) / num_processes ) * powers.size()) );

        size_t stone_block_start = size_t( floor( ( floor( pid   / double(processors_per_likelihood)) / (double(num_processes) / processors_per_likelihood) ) * powers.size() ) );
        size_t stone_block_end   = size_t( floor( ( ceil( (pid + 1) / double(processors_per_likelihood)) / (double(num_processes) / processors_per_likelihood) ) * powers.size() ) );
        
        for (size_t i = stone_block_start; i < stone_block_end; ++i)
        {
            runStone( i, gen, burnin_fraction, pre_burnin_generations, tuning_interval, false, checkpoint_file, checkpoint_interval );
        }
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


void PowerPosteriorAnalysis::runStone(size_t idx, size_t gen, double burnin_fraction, size_t pre_burnin_generations, size_t tuning_interval, bool one_only, const path &checkpoint_file, size_t checkpoint_interval)
{
    // create the directory if necessary
    if (filename.filename().empty() or filename.filename() == "." or filename.filename() == "..")
    {
        throw RbException("Please provide a filename with an extension");
    }

    std::string stone_tag = "_stone_" + std::to_string(idx + 1); // number the stones from 1 to n, not from 0 to (n - 1)

    path stoneFileName = appendToStem(filename, stone_tag);
    createDirectoryForFile( stoneFileName );

    // reset the sampler
    sampler->reset();

    size_t preburninPrintInterval = size_t( round( fmax(1, pre_burnin_generations/40.0) ) );
    size_t printInterval = size_t( round( fmax(1, gen/40.0) ) );
    
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
        
        if (num_processes == 1 or one_only == true or resume_from_checkpoint == true)
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
    
    sampler->addFileMonitorExtension(stone_tag, false);
    
    size_t k_start = 0;
    size_t burnin = 0;
    
    if ( resume_from_checkpoint )
    {
        // If we call .runOneStone() after .initializeFromCheckpoint(), we need the user to provide an index for which we do in fact
        // have a checkpoint file (whose name is now stored in ckp_stone_file). Trying to call the method with a different index
        // will throw an error.
        auto it = ckp_stone_file.find( idx );
        if ( it == ckp_stone_file.end() )
        {
            throw RbException() << "Stone " << (idx + 1) << " was not initialized from checkpoint.";
        }
        
        // We checked that the checkpoint file exists in initializeFromCheckpoint(); now we load it. We do not check again here:
        // this amounts to the assumption that it has not been deleted between the .initializeFromCheckpoint() and .runStone() calls.
        const path &stone_file_name = it->second;

        // give the stone back its correct power and planned_burnin
        double stone_power;
            
        // assemble the new filename
        path stone_info_file_name = appendToStem(stone_file_name, "_stone_info");

        // open file and initialize variables for parsing
        std::ifstream in_file_stone_info( stone_info_file_name.string() );
        std::string line;
        std::map<std::string, std::string> pars;
        
        // command-processing loop
        while ( in_file_stone_info.good() )
        {
            // read a line
            safeGetline( in_file_stone_info, line );
            
            if ( line != "" )
            {
                std::vector<std::string> key_value;
                StringUtilities::stringSplit(line, " = ", key_value);
                pars.insert( std::pair<std::string, std::string>(key_value[0], key_value[1]) );
            }
            
        }
        
        stone_power = StringUtilities::asDoubleNumber( pars["power"] );
        burnin = size_t( StringUtilities::asIntegerNumber( pars["planned_burnin"] ) );
        
        // clean up
        in_file_stone_info.close();

        sampler->setLikelihoodHeat( stone_power );

        // We only restore from checkpoint *after* we have set the correct heat, so that the DAG touches use the right likelihood
        // temperature. This is different from setting the chain posterior heat in MC^3, which only takes effect when applying moves.
        sampler->setCheckpointFile( stone_file_name );
        sampler->initializeSamplerFromCheckpoint();
        
        k_start = sampler->getCurrentGeneration();

        if (pid == pid_to_print and k_start < burnin)
        {
            RBOUT( "   Stone " + std::to_string( idx + 1 ) + ": finishing a previously scheduled and interrupted burnin stage; burninFraction of the current run() is ignored." );
        }
    }
    else
    {
        sampler->setLikelihoodHeat( powers[idx] );
        
        // Let's do a stone-specific pre-burnin: this is different from the pre-burnin performed by .burnin(), which is global.
        for (size_t k=1; k<=pre_burnin_generations; k++)
        {
            if (pid == pid_to_print)
            {
                if (k % preburninPrintInterval == 0)
                {
                    std::cout << "-";
                    std::cout.flush();
                }
            }

            sampler->nextCycle(false);
            
            // check for autotuning
            if ( k % tuning_interval == 0 && k != pre_burnin_generations )
            {
                sampler->tune();
            }
        }

        burnin = size_t( ceil( burnin_fraction * gen ) );
    }

    // We will run this stone for an additional 'gen' iterations, even if we are resuming from checkpoint.
    const size_t k_final = k_start + gen;
    
    // Monitor. Note that if we are running just one stone at a time, only one process is allowed to write the monitors.
    if (not one_only or (one_only and pid == pid_to_print))
    {
        sampler->startMonitors(gen, resume_from_checkpoint);
        if ( not resume_from_checkpoint )
        {
            sampler->writeMonitorHeaders( false );
        }
        // Do not re-print the checkpoint generation before the next cycle (avoids duplicate monitor lines).
        if ( not resume_from_checkpoint )
        {
            sampler->monitor(k_start);
        }
    }
    
    // Open the output file: append if resuming, overwrite if starting fresh.
    std::ofstream outStream;
    if ( resume_from_checkpoint )
    {
        outStream.open( stoneFileName.string(), std::ios::app );
    }
    else
    {
        outStream.open( stoneFileName.string() );
        outStream << "state\t" << "power\t" << "likelihood" << std::endl;
    }
    
    double p = powers[idx];

    // We don't want the sampler's internal generation counter to keep accumulating across stones: it should stay in sync
    // with the stone-local counter k. This is especially important for checkpointing. On a fresh start k_start is 0;
    // on resumption it holds the checkpointed iteration. Then nextCycle(true) increments generation each call,
    // keeping it equal to k.
    sampler->setCurrentGeneration( k_start );

    for (size_t k = k_start + 1; k <= k_final; ++k)
    {
        
        if (pid == pid_to_print)
        {
            if (k % printInterval == 0)
            {
                if (k <= burnin)
                {
                    std::cout << "+";
                    std::cout.flush();
                }
                else
                {
                    std::cout << "*";
                    std::cout.flush();
                }
            }
        }
        
        sampler->nextCycle( true );

        // monitor
        if (not one_only or (one_only and pid == pid_to_print))
        {
            sampler->monitor(k);
        }
        
        // sample the likelihood
        if ( k > burnin && k % sampleFreq == 0 )
        {
            // compute the joint likelihood
            double likelihood = sampler->getModelLnProbability(true);
            outStream << k << "\t" << p << "\t" << likelihood << std::endl;
        }
        
        // periodically checkpoint (including during burnin)
        if ( checkpoint_interval != 0 && (k % checkpoint_interval) == 0 )
        {
            checkpoint( idx, checkpoint_file, burnin );
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
        std::string stone_tag = "_stone_" + std::to_string(idx + 1); // number the stones from 1 to n, not from 0 to (n - 1)
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

