#ifndef AlleleFrequencySimulator_H
#define AlleleFrequencySimulator_H

#include <vector>

#include "HomologousDiscreteCharacterData.h"
#include "PoMoState.h"

namespace RevBayesCore {
class Tree;
    
    /**
     * This class provides an allele frequency simulator.
     *
     * We simulate the allele frequencies under a Moran process along a phylogeny.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2022-05-17, version 1.1
     */
    class AlleleFrequencySimulator {
        
    public:
        
        AlleleFrequencySimulator(Tree* t, const std::vector<long>& ps, double gt, size_t ns, const std::vector<double>& mr, const std::vector<long>& s, double r);
        
        void                                            simulateAlleleFrequencies( const std::string& fn, bool only_Variable ) const;
        
    private:
        
        bool                                            isVariable(std::vector<int>& site_pattern) const;
//        void                                            simulate( const TopologyNode& n, const std::vector<long>& states, std::vector<std::vector<int> >& taxa, std::vector<bool>& monomorphic ) const;
        bool                                            simulate(const TopologyNode& n, long states, std::vector<int>& site_pattern, bool& monomorphic ) const;
        long                                            simulateAlongBranch(size_t root_index, long root_start_state, double branch_length ) const;
//        std::vector<long>                               simulateAlongBranch( size_t root_index, const std::vector<long>& root_start_states, double branch_length ) const;
        void                                            writeCountsFile(const std::string& fn, const std::vector<std::vector<int> >& taxa) const;
        
        Tree*                                           tree;
        std::vector<long>                               population_sizes;
        double                                          generation_time;
        size_t                                          num_sites;
        std::vector<double>                             mutation_rates;
        std::vector<long>                               samples_per_species;
        double                                          root_branch;
        
        
    };
    
}


#endif
