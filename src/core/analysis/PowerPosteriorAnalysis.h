#ifndef PowerPosteriorAnalysis_H
#define PowerPosteriorAnalysis_H

#include "Cloneable.h"
#include "Parallelizable.h"
#include "RbVector.h"
#include "RbFileManager.h"

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
    class PowerPosteriorAnalysis : public Cloneable, public Parallelizable {
        
    public:
        PowerPosteriorAnalysis(MonteCarloSampler *m, const path &fn, size_t k);
        PowerPosteriorAnalysis(const PowerPosteriorAnalysis &a);
        virtual                                ~PowerPosteriorAnalysis(void);                               //!< Virtual destructor
        
        PowerPosteriorAnalysis&                 operator=(const PowerPosteriorAnalysis &a);
        
        // public methods
        PowerPosteriorAnalysis*                 clone(void) const;
        void                                    burnin(size_t g, size_t ti);
        void                                    checkpoint(size_t stone_idx, const path &bcp_file, size_t planned_burnin) const;
        std::vector<double>                     getPowers(void) const;
        void                                    initializeFromCheckpoint(const path &bcp_file, const std::vector<size_t> &stone_indices); //!< Makes sure checkpoint files exist but does not load the sampler
        void                                    runAll(size_t g, double burn_frac, size_t preburn_gen, size_t tune_int, const path &cp_file, size_t ci = 0);
        void                                    runStone(size_t idx, size_t g, double burn_frac, size_t preburn_gen, size_t tune_int, bool one_only, const path &cp_file, size_t ci = 0);
        void                                    summarizeStones(void);
        void                                    setPowers(const std::vector<double> &p);
        void                                    setSampleFreq(size_t sf);
        
    private:
        
        void                                    initMPI(void);
        
        // members
        path                                    filename;
        std::vector<double>                     powers;
        MonteCarloSampler*                      sampler;
        size_t                                  sampleFreq;
        size_t                                  processors_per_likelihood;
        
        bool                                    resume_from_checkpoint;                                     //!< Set true by initializeFromCheckpoint()
        std::map<size_t, path>                  ckp_stone_file;                                             //!< Ordered by increasing stone index.
    };
    
}

#endif
