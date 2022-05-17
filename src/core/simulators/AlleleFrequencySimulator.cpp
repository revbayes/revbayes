#include "AlleleFrequencySimulator.h"
#include "DiscreteTaxonData.h"
#include "DistributionBinomial.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "Tree.h"

using namespace RevBayesCore;

AlleleFrequencySimulator::AlleleFrequencySimulator(Tree* t, const std::vector<long>& ps, double gt, size_t ns, const std::vector<double>& mr, const std::vector<long>& s, double r) :
tree( t ),
population_sizes( ps ),
generation_time( gt ),
num_sites( ns ),
mutation_rates( mr ),
samples_per_species( s ),
root_branch( r )
{
    
}


void AlleleFrequencySimulator::simulateAlleleFrequencies( const std::string& fn ) const
{
    // first, get some variables/settings for the simulation
    size_t num_tips     = tree->getNumberOfTips();
    size_t root_index   = tree->getRoot().getIndex();
    
    RandomNumberGenerator* rng = GLOBAL_RNG;

    // create a vector of taxon data
    std::vector<std::vector<int> > taxa = std::vector<std::vector<int> >( num_tips, std::vector<int>() );
    std::vector<long> root_start_states =  std::vector<long>( num_sites, 0);
    for ( size_t i = 0; i < num_sites; ++i )
    {

        // draw the state
        double u = rng->uniform01();
        if ( u < mutation_rates[0] / ( mutation_rates[0] + mutation_rates[1] ) )
        {
            root_start_states[i] = 0;
        }
        else
        {
            root_start_states[i] = population_sizes[root_index];
        }
        
    }
    // simulate the root sequence
    std::vector<long> root_states = simulateAlongBranch( root_index, root_start_states, root_branch );

    // recursively simulate the sequences
    simulate( tree->getRoot(), root_states, taxa );
    
    writeCountsFile( fn, taxa );
    
}


void AlleleFrequencySimulator::simulate( const TopologyNode& n, const std::vector<long>& states, std::vector<std::vector<int> >& taxa ) const
{
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    if ( n.isTip() == true )
    {
        //
        size_t node_index = n.getIndex();
        double this_populuation_size = population_sizes[node_index];
        double samples = samples_per_species[node_index];
        
        std::vector<int> &taxon = taxa[node_index];
        
        for ( size_t i = 0; i < num_sites; ++i )
        {
            int tip_sample = RbStatistics::Binomial::rv(samples, states[i]/this_populuation_size, *rng);
            taxon.push_back( tip_sample );
        }
    }
    else
    {
        const TopologyNode& left = n.getChild(0);
        size_t left_index  = left.getIndex();
        double left_branch = left.getBranchLength();
        std::vector<long> left_states = simulateAlongBranch( left_index, states, left_branch );
        simulate( left, left_states, taxa );
        
        // the right child
        const TopologyNode& right = n.getChild(1);
        size_t right_index  = right.getIndex();
        double right_branch = right.getBranchLength();
        std::vector<long> right_states = simulateAlongBranch( right_index, states, right_branch );
        simulate( right, right_states, taxa );
    }
}



std::vector<long> AlleleFrequencySimulator::simulateAlongBranch( size_t node_index, const std::vector<long>& root_start_states, double branch_length ) const
{
    
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    double current_time = 0.0;
    double this_populuation_size = population_sizes[node_index];
    double per_generation_mutation_rate_0 = mutation_rates[0] * this_populuation_size / generation_time;
    double per_generation_mutation_rate_1 = mutation_rates[1] * this_populuation_size / generation_time;
    
    std::vector<long> current_states = root_start_states;
    
    while ( current_time < branch_length )
    {
        current_time += generation_time;
        
        for ( size_t i = 0; i < num_sites; ++i )
        {
            // only perform drift if we are in a polymorphic state
            if ( current_states[i] > 0 && current_states[i] < this_populuation_size )
            {
                current_states[i] = RbStatistics::Binomial::rv(this_populuation_size, current_states[i]/this_populuation_size, *rng);
            }
            
            // check for the boundary states at which mutation happen
            if ( current_states[i] == 0 )
            {
                int num_mutants = RbStatistics::Binomial::rv(this_populuation_size, per_generation_mutation_rate_0, *rng);
                current_states[i] += num_mutants;
            }
            else if ( current_states[i] == this_populuation_size )
            {
                int num_mutants = RbStatistics::Binomial::rv(this_populuation_size, per_generation_mutation_rate_1, *rng);
                current_states[i] -= num_mutants;
            }
        }
    }
    
    return current_states;
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
            int sampled_ones = taxa[tip_index][site];
            int samples = int(samples_per_species[tip_index]);
            out_stream << " " << (samples-sampled_ones) << "," << sampled_ones;
            
        }
        out_stream << std::endl;
        
    }
    
    // close the stream
    out_stream.close();
}

