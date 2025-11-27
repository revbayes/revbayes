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
     * We simulate the allele frequencies under a Moran process astd::int64_t a phylogeny.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2022-05-17, version 1.1
     */
    class AlleleFrequencySimulator : public Parallelizable, public Cloneable {
        
    public:
        
        AlleleFrequencySimulator(double gt, const std::vector<double>& mr, bool mg);
        
        AlleleFrequencySimulator*                       clone(void) const;
        
        void                                            simulateAlleleFrequencies(const Tree* t, const std::vector<std::int64_t>& ps, size_t ns, const std::vector<std::int64_t>& s, double r, const std::string& fn, bool only_Variable ) const;
        MatrixReal*                                     simulateAlleleFrequenciesMatrix( double t, std::int64_t ps, std::int64_t r ) const;
        RbVector<double>*                               simulateAlleleFrequenciesVector( double t, std::int64_t ps, std::int64_t r, size_t s ) const;
        RbVector<double>*                               simulateAlleleFrequenciesVectorEpoch( const std::vector<double>& t, const std::vector<std::int64_t>& ps, std::int64_t r, size_t s, std::int64_t f ) const;

    private:
        
        bool                                            isVariable(std::vector<int>& site_pattern) const;
        bool                                            simulateAlignment(const TopologyNode& n, std::int64_t states, const std::vector<std::int64_t>& ps, const std::vector<std::int64_t>& s, std::vector<int>& site_pattern, bool& monomorphic ) const;
        std::int64_t                                    simulateAlongBranch(double ps, std::int64_t root_start_state, double branch_length ) const;
        void                                            writeCountsFile(const Tree* t, const path& fn, const std::vector<std::vector<int> >& taxa, const std::vector<std::int64_t>& s) const;
        
        double                                          generation_time;
        std::vector<double>                             mutation_rates;
        bool                                            moran_generations;
        
        
    };
    
}


#endif
