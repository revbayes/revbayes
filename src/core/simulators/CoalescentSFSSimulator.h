#ifndef CoalescentSFSSimulator_H
#define CoalescentSFSSimulator_H

#include <vector>

#include "HomologousDiscreteCharacterData.h"
#include "Parallelizable.h"
#include "PoMoState.h"

namespace RevBayesCore {
class RandomNumberGenerator;
    
    /**
     * This class provides an allele frequency simulator.
     *
     * We simulate the allele frequencies under a Moran process along a phylogeny.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2022-05-17, version 1.1
     */
    class CoalescentSFSSimulator : public Parallelizable, public Cloneable {
        
    public:
        
        CoalescentSFSSimulator(const std::vector<double>& p, const std::vector<double>& cp, double gt, double mr);
        
        CoalescentSFSSimulator*                         clone(void) const;
        
        RbVector<long>*                                 simulateSFS(long s, long r ) const;

    private:
        
        double                                          simulateCoalescentTime(double t, size_t a, RandomNumberGenerator* r) const;
        size_t                                          simulateMutations(size_t i, long s, const std::vector<std::vector<size_t> >& c, const std::vector<double>& a, std::vector<long>& t, RandomNumberGenerator* r) const;

        double                                          generation_time;
        double                                          mutation_rate;
        std::vector<double>                             population_sizes;
        std::vector<double>                             change_points;
        
        
    };
    
}


#endif
