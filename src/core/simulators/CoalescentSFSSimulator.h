#ifndef CoalescentSFSSimulator_H
#define CoalescentSFSSimulator_H

#include <vector>

#include "DemographicFunction.h"
#include "Parallelizable.h"
#include "RbVector.h"

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
        
        CoalescentSFSSimulator(const RbVector<DemographicFunction>& p, const std::vector<double>& cp, double gt, const std::string& pl);
        
        CoalescentSFSSimulator*                         clone(void) const;
        
        RbVector<long>*                                 simulateSFS(double mr, long s, long r) const;
        void                                            simulateCoalescent(long s, long r, const path& f_stats, const path& f_trees) const;

    private:
        
        double                                          simulateCoalescentTime(double t, size_t a, RandomNumberGenerator* r) const;
        size_t                                          simulateMutations(double mr, size_t i, long s, const std::vector<std::vector<size_t> >& c, const std::vector<double>& a, std::vector<long>& t, RandomNumberGenerator* r) const;

        double                                          generation_time;
        RbVector<DemographicFunction>                   demographies;
        std::vector<double>                             change_points;
        double                                          ploidy_factor;
        
        
    };
    
}


#endif
