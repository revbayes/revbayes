#ifndef AlleleFrequencySimulator_H
#define AlleleFrequencySimulator_H

#include <vector>

#include "HomologousDiscreteCharacterData.h"
#include "Parallelizable.h"
#include "PoMoState.h"

namespace RevBayesCore {
class Tree;
class TransitionProbabilityMatrix;
class MatrixReal;
    
    /**
     * This class provides an allele frequency simulator.
     *
     * We simulate the allele frequencies under a Moran process along a phylogeny.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2022-05-17, version 1.1
     */
    class AlleleFrequencySimulator : public Parallelizable, public Cloneable {
        
    public:
        
        AlleleFrequencySimulator(double gt, const std::vector<double>& mr, bool mg);
        
        AlleleFrequencySimulator*                       clone(void) const;
        
        void                                            simulateAlleleFrequencies(const Tree* t, const std::vector<long>& ps, size_t ns, const std::vector<long>& s, double r, const std::string& fn, bool only_Variable ) const;
        MatrixReal*                                     simulateAlleleFrequenciesMatrix( double t, long ps, long r ) const;
        RbVector<double>*                               simulateAlleleFrequenciesVector( double t, long ps, long r, size_t s ) const;

    private:
        
        bool                                            isVariable(std::vector<int>& site_pattern) const;
        bool                                            simulateAlignment(const TopologyNode& n, long states, const std::vector<long>& ps, const std::vector<long>& s, std::vector<int>& site_pattern, bool& monomorphic ) const;
        long                                            simulateAlongBranch(double ps, long root_start_state, double branch_length ) const;
        void                                            writeCountsFile(const Tree* t, const std::string& fn, const std::vector<std::vector<int> >& taxa, const std::vector<long>& s) const;
        
        double                                          generation_time;
        std::vector<double>                             mutation_rates;
        bool                                            moran_generations;
        
        
    };
    
}


#endif
