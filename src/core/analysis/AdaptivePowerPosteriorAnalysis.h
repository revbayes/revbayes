#ifndef AdaptivePowerPosteriorAnalysis_H
#define AdaptivePowerPosteriorAnalysis_H

#include "Cloneable.h"
#include "Parallelizable.h"
#include "RbVector.h"

namespace RevBayesCore {
    
    class MonteCarloSampler;
    
    /**
     * @brief Power posterior analysis class.
     *
     * A power posterior analysis runs an analysis for a vector of powers
     * where the likelihood during each analysis run is raised to the given power.
     * The likelihood values and the current powers are stored in a file.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2012-06-17
     *
     */
    class AdaptivePowerPosteriorAnalysis : public Cloneable, public Parallelizable {
        
    public:
        AdaptivePowerPosteriorAnalysis(MonteCarloSampler *m, const std::string &fn, size_t k);
        AdaptivePowerPosteriorAnalysis(const AdaptivePowerPosteriorAnalysis &a);
        virtual                                ~AdaptivePowerPosteriorAnalysis(void);                               //!< Virtual destructor
        
        AdaptivePowerPosteriorAnalysis&         operator=(const AdaptivePowerPosteriorAnalysis &a);
        
        // public methods
        AdaptivePowerPosteriorAnalysis*         clone(void) const;
        void                                    burnin(size_t g, size_t ti);
        void                                    runAll(size_t g, size_t bg, size_t ti);
        void                                    summarizeStones(void);
        void                                    setSampleFreq(size_t sf);
        void                                    setTolerance(double tol);
        
    private:
        
        void                                    initMPI(void);
        
        void                                    runStone(size_t idx, size_t g, size_t bg, size_t ti);
        double                                  adapt(void);
        void                                    addStones(void);
        
        double                                  trapezoidRule(double fa, double fb, double a, double b);
        double                                  simpsonsRule(double fa, double fb, double fc, double a, double b, double c);
        
        // members
        std::string                             filename;
        
        std::vector<double>                     powers;
        std::map<double, double>                average_likelihoods;
        std::map<double, bool>                  dirty_stones;
        std::map<double, size_t>                indices;
        double                                  tolerance;
        double                                  old_marginal_likelihood;
        
        size_t                                  worst_stone;
        double                                  new_beta;
        
        MonteCarloSampler*                      sampler;
        size_t                                  sampleFreq;                                                                     //!< The rate of the distribution
        size_t                                  processors_per_likelihood;

    };
    
}

#endif
