#include "AlleleFrequencySimulator.h"
#include "DiscreteTaxonData.h"
#include "DistributionBinomial.h"
#include "MatrixReal.h"
#include "ProgressBar.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "TransitionProbabilityMatrix.h"
#include "Tree.h"

using namespace RevBayesCore;

AlleleFrequencySimulator::AlleleFrequencySimulator(Tree* t, const std::vector<long>& ps, double gt, size_t ns, const std::vector<double>& mr, const std::vector<long>& s, double r, bool mg) :
tree( t ),
population_sizes( ps ),
generation_time( gt ),
num_sites( ns ),
mutation_rates( mr ),
samples_per_species( s ),
root_branch( r ),
moran_generations( mg )
{
    
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



void AlleleFrequencySimulator::simulateAlleleFrequencies( const std::string& fn, bool only_variable ) const
{
    // first, get some variables/settings for the simulation
    size_t num_tips     = tree->getNumberOfTips();
    size_t root_index   = tree->getRoot().getIndex();
    
    RandomNumberGenerator* rng = GLOBAL_RNG;

    // create a vector of taxon data
    std::vector<std::vector<int> > taxa = std::vector<std::vector<int> >( num_sites, std::vector<int>() );
    std::vector<bool> monomorphic =  std::vector<bool>( num_sites, true);
    
    // start the progress bar
    ProgressBar progress = ProgressBar(num_sites, 0);
    
    // Print progress bar (68 characters wide)
    progress.start();
    
    size_t num_attempts = 0;
    
    std::vector<long> root_start_states =  std::vector<long>( num_sites, 0);
    for ( size_t i = 0; i < num_sites; ++i )
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
            success = simulate( tree->getRoot(), root_state, site_pattern, mono );
            success = only_variable == false || isVariable( site_pattern );
            
            
        } while ( success == false );
        
        taxa[i] = site_pattern;
        monomorphic[i] = mono;
        
        progress.update(i);
        
    }
    
    progress.finish();

    
    size_t num_monomorphic = 0;
    for ( size_t i = 0; i < num_sites; ++i )
    {
        if ( monomorphic[i] == false )
        {
            ++num_monomorphic;
        }
    }
    std::cerr << "#Monomorphic sites:\t\t" << (num_sites-num_monomorphic) << std::endl;
    std::cerr << "#Biallelic sites:\t\t" << num_monomorphic << std::endl;
    std::cerr << "#Invariant sites:\t\t" << (num_attempts-num_sites) << std::endl;
    
    
    writeCountsFile( fn, taxa );
    
}


MatrixReal* AlleleFrequencySimulator::simulateAlleleFrequenciesMatrix( double time, long population_size, long reps ) const
{
    
    MatrixReal* tpm = new MatrixReal(population_size+1);
    
    // start the progress bar
    ProgressBar progress = ProgressBar(population_size, 0);
    
    // Print progress bar (68 characters wide)
    progress.start();
    
    for ( long i=0; i<=population_size; ++i )
    {
        
        std::vector<long> counts = std::vector<long>(population_size+1, 0);
        
        for (size_t r=0; r<reps; ++r)
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
    
    return tpm;
}



bool AlleleFrequencySimulator::simulate( const TopologyNode& n, long state, std::vector<int>& taxa, bool& monomorphic ) const
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
        simulate( left, left_state, taxa, monomorphic );
        
        // the right child
        const TopologyNode& right = n.getChild(1);
        size_t right_index  = right.getIndex();
        double right_branch = right.getBranchLength();
        long right_state = simulateAlongBranch( population_sizes[right_index], state, right_branch );
        simulate( right, right_state, taxa, monomorphic );
    }
    
    return true;
}



long AlleleFrequencySimulator::simulateAlongBranch( double this_populuation_size, long root_start_state, double branch_length ) const
{
    
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    double current_time = 0.0;
    double this_populuation_size = population_sizes[node_index];
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
                
            double u = rng->uniform01();
            if ( per_generation_mutation_rate_0 > u )
            {
                current_state = 1;
            }
        }
        else if ( current_state == this_populuation_size )
        {
//            int num_mutants = RbStatistics::Binomial::rv(this_populuation_size, per_generation_mutation_rate_1, *rng);
//            current_state -= num_mutants;
                
            double u = rng->uniform01();
            if ( per_generation_mutation_rate_1 > u )
            {
                    current_state = this_populuation_size-1;
            }
        }
    }
    
    return current_state;
}


void AlleleFrequencySimulator::writeCountsFile(const std::string& fn, const std::vector<std::vector<int> >& taxa) const
{
    
    // first, get some variables/settings for the simulation
    size_t num_tips     = tree->getNumberOfTips();
    
    // the filestream object
    std::fstream out_stream;

    RbFileManager f = RbFileManager(fn);
    f.createDirectoryForFile();

    // open the stream to the file
    out_stream.open( f.getFullFileName().c_str(), std::fstream::out );
    
    /*
     COUNTSFILE NPOP 12 NSITES 1000
     CHROM POS Gorilla_beringei Gorilla_gorilla ...
     chr1 41275799 6,0,0,0 2,0,0,0 ...
     */
    out_stream << "COUNTSFILE NPOP " << num_tips << " NSITES " << num_sites << std::endl;
    out_stream << "CHROM POS";
    for (size_t i=0; i<num_tips; ++i)
    {
        out_stream << " " << tree->getTipNode(i).getName();
    }
    out_stream << std::endl;
    for (size_t site=0; site<num_sites; ++site)
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
}

