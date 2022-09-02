#include "AlleleFrequencySimulator.h"
#include "DiscreteTaxonData.h"
#include "DistributionBinomial.h"
#include "DistributionGeometric.h"
#include "MatrixReal.h"
#include "MpiUtilities.h"
#include "ProgressBar.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "TransitionProbabilityMatrix.h"
#include "Tree.h"

#ifdef RB_MPI
#include <mpi.h>
#include <thread>
#include <chrono>
#endif

using namespace RevBayesCore;

AlleleFrequencySimulator::AlleleFrequencySimulator(double gt, const std::vector<double>& mr, bool mg) :
generation_time( gt ),
mutation_rates( mr ),
moran_generations( mg )
{
    
}


AlleleFrequencySimulator* AlleleFrequencySimulator::clone( void ) const
{
    return new AlleleFrequencySimulator( *this );
}


bool AlleleFrequencySimulator::isVariable(std::vector<int>& site_pattern) const
{
    
    int ref = site_pattern[0];
    for (size_t i=1; i<site_pattern.size(); ++i)
    {
        if ( ref != site_pattern[i] )
        {
            return true;
        }
    }
    
    return false;    
}



void AlleleFrequencySimulator::simulateAlleleFrequencies( const Tree* tree, const std::vector<long>& population_sizes, size_t num_sites, const std::vector<long>& samples_per_species, double root_branch, const std::string& fn, bool only_variable ) const
{
    // first, get some variables/settings for the simulation
    size_t num_tips     = tree->getNumberOfTips();
    size_t root_index   = tree->getRoot().getIndex();
    
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    // forward the rng for different processes
#ifdef RB_MPI
    for ( size_t i=active_PID; i<pid; ++i)
    {
        // we fast forward 7 times, just to be sure
        rng->uniform01();
        rng->uniform01();
        rng->uniform01();
        rng->uniform01();
        rng->uniform01();
        rng->uniform01();
        rng->uniform01();
    }
#endif
    
    size_t num_sites_this_process = num_sites;
#ifdef RB_MPI
    num_sites_this_process = num_sites / num_processes;
#endif
    

    // create a vector of taxon data
    std::vector<std::vector<int> > taxa = std::vector<std::vector<int> >( num_sites_this_process, std::vector<int>() );
    std::vector<bool> monomorphic =  std::vector<bool>( num_sites_this_process, true);
    
    // start the progress bar
    ProgressBar progress = ProgressBar(num_sites_this_process, 0);
    
    // Print progress bar (68 characters wide)
    progress.start();
    
    size_t num_attempts = 0;
    
    for ( size_t i = 0; i < num_sites_this_process; ++i )
    {
        std::vector<int> site_pattern = std::vector<int>( num_tips, 0);
        bool mono = true;

        bool success = false;
        do {
            ++num_attempts;
            long root_start_state = 0;
            
            // draw the state
            double u = rng->uniform01();
            if ( u < mutation_rates[0] / ( mutation_rates[0] + mutation_rates[1] ) )
            {
                root_start_state = 0;
            }
            else
            {
                root_start_state = population_sizes[root_index];
            }
            
            // simulate the root sequence
            long root_state = simulateAlongBranch( population_sizes[root_index], root_start_state, root_branch );

            // recursively simulate the sequences
            mono = true;
            success = simulateAlignment( tree->getRoot(), root_state, population_sizes, samples_per_species, site_pattern, mono );
            success = only_variable == false || isVariable( site_pattern );
            
            
        } while ( success == false );
        
        taxa[i] = site_pattern;
        monomorphic[i] = mono;
        
        progress.update(i);
        
    }
    
    progress.finish();

    
    size_t num_monomorphic = 0;
    for ( size_t i = 0; i < num_sites_this_process; ++i )
    {
        if ( monomorphic[i] == false )
        {
            ++num_monomorphic;
        }
    }
#ifdef RB_MPI
    MPI_Barrier( MPI_COMM_WORLD );

    // we only need to send message if there is more than one process
    if ( num_processes > 1 )
    {

        // send the number of monomorphic sites from the helpers to the master
        if ( process_active == false )
        {
            // send from the workers the number of monomorphic sites to the master
            MPI_Send(&num_monomorphic, 1, MPI_INT, active_PID, 0, MPI_COMM_WORLD);
        }

        // receive the number of monomorphic sites from the helpers
        if ( process_active == true )
        {
            for (size_t i=active_PID+1; i<active_PID+num_processes; ++i)
            {
                int tmp = 0;
                MPI_Status status;
                MPI_Recv(&tmp, 1, MPI_INT, int(i), 0, MPI_COMM_WORLD, &status);
                num_monomorphic += tmp;
            }
        }
    }
#endif
    
    
#ifdef RB_MPI
    // we only need to send message if there is more than one process
    if ( process_active == true )
    {
#endif
    std::cerr << "#Monomorphic sites:\t\t" << (num_sites-num_monomorphic) << std::endl;
    std::cerr << "#Biallelic sites:\t\t" << num_monomorphic << std::endl;
#ifdef RB_MPI
    }
#endif
    
    writeCountsFile( tree, fn, taxa, samples_per_species );
    
    
#ifdef RB_MPI
    MpiUtilities::synchronizeRNG( MPI_COMM_WORLD );
#endif
    
}


MatrixReal* AlleleFrequencySimulator::simulateAlleleFrequenciesMatrix( double time, long population_size, long reps ) const
{
    
    MatrixReal* tpm = new MatrixReal(population_size+1);
    
    // forward the rng for different processes
#ifdef RB_MPI
    RandomNumberGenerator* rng = GLOBAL_RNG;

    for ( size_t i=active_PID; i<pid; ++i)
    {
        // we fast forward 7 times, just to be sure
        rng->uniform01();
        rng->uniform01();
        rng->uniform01();
        rng->uniform01();
        rng->uniform01();
        rng->uniform01();
        rng->uniform01();
    }
#endif
    
    size_t reps_this_process = reps;
#ifdef RB_MPI
    reps_this_process = reps / num_processes;
#endif
    
    // start the progress bar
    ProgressBar progress = ProgressBar(population_size, 0);
    
    // Print progress bar (68 characters wide)
    progress.start();
    
    for ( long i=0; i<=population_size; ++i )
    {
        
        std::vector<long> counts = std::vector<long>(population_size+1, 0);
        
        for (size_t r=0; r<reps_this_process; ++r)
        {
            long obs_state = simulateAlongBranch(population_size, i, time);
            ++counts[obs_state];
        }
        
        for ( size_t s=0; s<=population_size; ++s )
        {
            (*tpm)[i][s] = double(counts[s])/reps;
        }
        
        progress.update(i);
        
    }
    
    progress.finish();
    
#ifdef RB_MPI
    MPI_Barrier( MPI_COMM_WORLD );
    
    // create a copy of the transition probability matrix
    MatrixReal tpm_backup = MatrixReal( *tpm );

    // we only need to send message if there is more than one process
    if ( num_processes > 1 )
    {

        std::vector< double > this_tpm_row = std::vector<double>(population_size+1,0.0);

        for (size_t i=active_PID; i<active_PID+num_processes; ++i)
        {
            
            for ( size_t s=0; s<=population_size; ++s )
            {
                if ( pid == i )
                {
                    this_tpm_row = tpm_backup[s];
                }
                
                MPI_Bcast(&this_tpm_row[0], population_size+1, MPI_DOUBLE, pid, MPI_COMM_WORLD);
                
                if ( pid != i )
                {
                    for ( size_t k=0; k<=population_size; ++k )
                    {
                        (*tpm)[s][k] += this_tpm_row[k];
                    }
                } // end-if non-sending process to add the transition probabilities
                
            } // end-for over all starting states for the transition probability matrix
            
        } // end-for over all processes
        
    } // end-if there are more than one process
#endif

    
    return tpm;
}


RbVector<double>* AlleleFrequencySimulator::simulateAlleleFrequenciesVector( double time, long population_size, long reps, size_t start_index ) const
{
    
    RbVector<double>* tpv = new RbVector<double>(population_size+1);
    
    // forward the rng for different processes
#ifdef RB_MPI
    RandomNumberGenerator* rng = GLOBAL_RNG;

    for ( size_t i=active_PID; i<pid; ++i)
    {
        // we fast forward 7 times, just to be sure
        rng->uniform01();
        rng->uniform01();
        rng->uniform01();
        rng->uniform01();
        rng->uniform01();
        rng->uniform01();
        rng->uniform01();
    }
#endif
    
    size_t reps_this_process = reps;
#ifdef RB_MPI
    reps_this_process = reps / num_processes;
#endif
    
    // start the progress bar
    ProgressBar progress = ProgressBar(reps_this_process, 0);
    
    // Print progress bar (68 characters wide)
    progress.start();
    
    std::vector<long> counts = std::vector<long>(population_size+1, 0);
        
    for (size_t r=0; r<reps_this_process; ++r)
    {
        long obs_state = simulateAlongBranch(population_size, start_index, time);
        ++counts[obs_state];
        
        progress.update(r);
    }
    progress.finish();
        
    for ( size_t s=0; s<=population_size; ++s )
    {
        (*tpv)[s] = double(counts[s])/reps;
    }
    
    
#ifdef RB_MPI
    MPI_Barrier( MPI_COMM_WORLD );
    
    // create a copy of the transition probability matrix
    RbVector<double> tpv_backup = RbVector<double>( *tpv );

    // we only need to send message if there is more than one process
    if ( num_processes > 1 )
    {
        std::vector< double > this_tpv_row = std::vector<double>(population_size+1,0.0);

        for (size_t i=active_PID; i<active_PID+num_processes; ++i)
        {
            
            if ( pid == i )
            {
                this_tpv_row = tpv_backup;
            }
                
            MPI_Bcast(&this_tpv_row[0], population_size+1, MPI_DOUBLE, pid, MPI_COMM_WORLD);
                
            if ( pid != i )
            {
                for ( size_t k=0; k<=population_size; ++k )
                {
                    (*tpv)[k] += this_tpv_row[k];
                }
            } // end-if non-sending process to add the transition probabilities
            
        } // end-for over all processes
        
    } // end-if there are more than one process
#endif

    
    return tpv;
}





RbVector<double>* AlleleFrequencySimulator::simulateAlleleFrequenciesVectorEpoch( const std::vector<double>& time, const std::vector<long>& population_size, long reps, size_t start_index, long final_population_size ) const
{
    // get some variables for the process
    
    // get the number of epochs
    size_t num_epochs = time.size();
    
    // create the vector with the counts for the final epoch
    RbVector<double>* tpv = new RbVector<double>(final_population_size+1);
    
    // get the gobal random number generator
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    // forward the rng for different processes
#ifdef RB_MPI
    for ( size_t i=active_PID; i<pid; ++i)
    {
        // we fast forward 7 times, just to be sure
        rng->uniform01();
        rng->uniform01();
        rng->uniform01();
        rng->uniform01();
        rng->uniform01();
        rng->uniform01();
        rng->uniform01();
    }
#endif
    
    size_t reps_this_process = reps;
#ifdef RB_MPI
    reps_this_process = reps / num_processes;
#endif
    
    // start the progress bar
    ProgressBar progress = ProgressBar(reps_this_process, 0);
    
    // Print progress bar (68 characters wide)
    progress.start();
    
    std::vector<long> counts = std::vector<long>(final_population_size+1, 0);
    
    // iterate over the replicates
    for (size_t r=0; r<reps_this_process; ++r)
    {
        
        // set the current state to the start state
        long current_state = start_index;
        
        // iterate over the epochs
        for (size_t epoch_index=0; epoch_index<num_epochs; ++epoch_index)
        {
            
            // get the population size for the current epoch
            long current_epoch_population_size = population_size[epoch_index];
            
            // if this is not the first epoch, then we need to rescale the index
            if ( epoch_index > 0 )
            {
                // get the population size for the current epoch
                long previous_epoch_population_size = population_size[epoch_index-1];
                double prev_freq = double(current_state)/previous_epoch_population_size;
                current_state = RbStatistics::Binomial::rv(current_epoch_population_size, prev_freq, *rng);
            }
            
            // get the population size for the current epoch
            double current_epoch_time = epoch_index > 0 ? time[epoch_index] - time[epoch_index-1] : time[epoch_index];
            
            // simulate the Moran process for the current epoch
            current_state = simulateAlongBranch(current_epoch_population_size, current_state, current_epoch_time);
        }
        
        // get the population size for the current epoch
        long previous_epoch_population_size = population_size[num_epochs-1];
        double prev_freq = double(current_state)/previous_epoch_population_size;
        long obs_state = RbStatistics::Binomial::rv(final_population_size, prev_freq, *rng);
        ++counts[obs_state];
        
        progress.update(r);
    }
    progress.finish();
        
    for ( size_t s=0; s<=final_population_size; ++s )
    {
        (*tpv)[s] = double(counts[s])/reps;
    }
    
    
#ifdef RB_MPI
    MPI_Barrier( MPI_COMM_WORLD );
    
    // create a copy of the transition probability matrix
    RbVector<double> tpv_backup = RbVector<double>( *tpv );

    // we only need to send message if there is more than one process
    if ( num_processes > 1 )
    {
        std::vector< double > this_tpv_row = std::vector<double>(final_population_size+1,0.0);

        for (size_t i=active_PID; i<active_PID+num_processes; ++i)
        {
            
            if ( pid == i )
            {
                this_tpv_row = tpv_backup;
            }
                
            MPI_Bcast(&this_tpv_row[0], final_population_size+1, MPI_DOUBLE, pid, MPI_COMM_WORLD);
                
            if ( pid != i )
            {
                for ( size_t k=0; k<=final_population_size; ++k )
                {
                    (*tpv)[k] += this_tpv_row[k];
                }
            } // end-if non-sending process to add the transition probabilities
            
        } // end-for over all processes
        
    } // end-if there are more than one process
#endif

    
    return tpv;
}



bool AlleleFrequencySimulator::simulateAlignment( const TopologyNode& n, long state, const std::vector<long>& population_sizes, const std::vector<long>& samples_per_species, std::vector<int>& taxa, bool& monomorphic ) const
{
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    if ( n.isTip() == true )
    {
        //
        size_t node_index = n.getIndex();
        double this_populuation_size = population_sizes[node_index];
        double samples = samples_per_species[node_index];
        
        int tip_sample = RbStatistics::Binomial::rv(samples, state/this_populuation_size, *rng);
        taxa[node_index] = tip_sample;
            
        if ( tip_sample > 0 && tip_sample < samples )
        {
            monomorphic = false;
        }
    }
    else
    {
        const TopologyNode& left = n.getChild(0);
        size_t left_index  = left.getIndex();
        double left_branch = left.getBranchLength();
        long left_state = simulateAlongBranch( population_sizes[left_index], state, left_branch );
        simulateAlignment( left, left_state, population_sizes, samples_per_species, taxa, monomorphic );
        
        // the right child
        const TopologyNode& right = n.getChild(1);
        size_t right_index  = right.getIndex();
        double right_branch = right.getBranchLength();
        long right_state = simulateAlongBranch( population_sizes[right_index], state, right_branch );
        simulateAlignment( right, right_state, population_sizes, samples_per_species, taxa, monomorphic );
    }
    
    return true;
}



long AlleleFrequencySimulator::simulateAlongBranch( double this_populuation_size, long root_start_state, double branch_length ) const
{
    
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    double current_time = 0.0;
    double this_generation_time = moran_generations ? generation_time / this_populuation_size : generation_time;
    double per_generation_mutation_rate_0 = mutation_rates[0] / generation_time;
    double per_generation_mutation_rate_1 = mutation_rates[1] / generation_time;
    
    long current_state = root_start_state;
    
    while ( current_time < branch_length )
    {
        current_time += this_generation_time;
        
        // only perform drift if we are in a polymorphic state
        if ( current_state > 0 && current_state < this_populuation_size )
        {
//            current_states = RbStatistics::Binomial::rv(this_populuation_size, current_states/this_populuation_size, *rng);
                
            double u = rng->uniform01();
            // pick a random ancestor individual
            if ( current_state/this_populuation_size > u )
            {
                // we picked an ancestor with a 1 state

                // pick a random individual to replace
                double v = rng->uniform01();
                if ( current_state/this_populuation_size < v )
                {
                    // we picked an individual with a 0 state, so we increase the counter of individuals with a 1
                    ++current_state;
                }
            }
            else
            {
                // we picked an ancestor with 0 state
                
                // pick a random individual to replace
                double v = rng->uniform01();
                if ( current_state/this_populuation_size > v )
                {
                    // we picked an individual with a 1 state, so we decrease the counter of individuals with a 1
                    --current_state;
                }
            }
        }
        // check for the boundary states at which mutation happen
        else if ( current_state == 0 )
        {
//            int num_mutants = RbStatistics::Binomial::rv(this_populuation_size, per_generation_mutation_rate_0, *rng);
//            current_state += num_mutants;
            
            // draw the waiting time until the next mutation
            int mutation_wait_generations = RbStatistics::Geometric::rv(per_generation_mutation_rate_0, *rng);
            
            // fast forward the generation time
            current_time += mutation_wait_generations * this_generation_time;
            
            // accept the mutation only if it happened within the time of this branch
            if ( current_time <= branch_length )
            {
                current_state = 1;
            }
        }
        else if ( current_state == this_populuation_size )
        {
//            int num_mutants = RbStatistics::Binomial::rv(this_populuation_size, per_generation_mutation_rate_1, *rng);
//            current_state -= num_mutants;
            
            // draw the waiting time until the next mutation
            int mutation_wait_generations = RbStatistics::Geometric::rv(per_generation_mutation_rate_1, *rng);
            
            // fast forward the generation time
            current_time += mutation_wait_generations * this_generation_time;
            
            // accept the mutation only if it happened within the time of this branch
            if ( current_time <= branch_length )
            {
                current_state = this_populuation_size-1;
            }
//            if ( per_generation_mutation_rate_1 > u )
//            {
//                current_state = this_populuation_size-1;
//            }
        }
    }
    
    return current_state;
}


void AlleleFrequencySimulator::writeCountsFile(const Tree* tree, const std::string& fn, const std::vector<std::vector<int> >& taxa, const std::vector<long>& samples_per_species) const
{
    
    // first, get some variables/settings for the simulation
    size_t num_tips     = tree->getNumberOfTips();
    size_t num_sites    = taxa.size();
#ifdef RB_MPI
    if ( num_processes > 1 )
    {
        num_sites *= num_processes;
    }
#endif
    
    // the filestream object
    std::fstream out_stream;

    RbFileManager f = RbFileManager(fn);
    
#ifdef RB_MPI
    if ( process_active == true )
    {
#endif
        f.createDirectoryForFile();
#ifdef RB_MPI
    }
#endif
    
    // open the stream to the file
    out_stream.open( f.getFullFileName().c_str(), std::fstream::out );
    
    /*
     COUNTSFILE NPOP 12 NSITES 1000
     CHROM POS Gorilla_beringei Gorilla_gorilla ...
     chr1 41275799 6,0,0,0 2,0,0,0 ...
     */
    
#ifdef RB_MPI
    MPI_Barrier( MPI_COMM_WORLD );
    if ( process_active == true )
    {
#endif
    out_stream << "COUNTSFILE NPOP " << num_tips << " NSITES " << num_sites << std::endl;
    out_stream << "CHROM POS";
    for (size_t i=0; i<num_tips; ++i)
    {
        out_stream << " " << tree->getTipNode(i).getName();
    }
    out_stream << std::endl;
#ifdef RB_MPI
    }
#endif
    
    size_t num_sites_per_process = taxa.size();
        
    
#ifdef RB_MPI
    for (size_t writing_pid=active_PID; writing_pid<active_PID+num_processes; ++writing_pid)
    {
        if ( writing_pid == pid )
        {
            // advance to end of
            out_stream.seekg(0, std::ios::end);
            
#endif
    for (size_t site=0; site<num_sites_per_process; ++site)
    {
        
        out_stream << "? ?";
        for (size_t i=0; i<num_tips; ++i)
        {
            size_t tip_index = tree->getTipNode(i).getIndex();
            int sampled_ones = taxa[site][tip_index];
            int samples = int(samples_per_species[tip_index]);
            out_stream << " " << (samples-sampled_ones) << "," << sampled_ones;
            
        }
        out_stream << std::endl;
        
    }
    
    // close the stream
    out_stream.close();
            
#ifdef RB_MPI
        }
        
        MPI_Barrier( MPI_COMM_WORLD );
    }
#endif
}

