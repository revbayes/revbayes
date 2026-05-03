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
        virtual                                ~PowerPosteriorAnalysis(void);                   //!< Virtual destructor
        
        PowerPosteriorAnalysis&                 operator=(const PowerPosteriorAnalysis &a);
        
        // public methods
        PowerPosteriorAnalysis*                 clone(void) const;
        void                                    burnin(size_t g, size_t ti, const path &cp_file, size_t ci = 0);
        void                                    checkpoint(const path &bcp_file) const;
        void                                    checkpoint(size_t stone_idx, const path &bcp_file, size_t planned_burnin) const;
        std::vector<double>                     getPowers(void) const;
        size_t                                  getStepNumber(void) const;                      //!< What is the largest number of stones to be executed by the same worker (i.e., a process or a group of processors_per_likelihood processes)?
        void                                    initializeFromCheckpoint(const path &bcp_file); //!< Makes sure checkpoint file exists but does not load the sampler
        void                                    initializeFromCheckpoint(const path &bcp_file, const std::vector<size_t> &stone_indices); //!< Does the same, but for each one of the specified stones
        void                                    initializeFromCheckpoint(const path &bcp_file, const std::vector<std::vector<size_t>> &stone_sequences_per_worker); //!< Records stone sequence per worker; only first stone per sequence needs a checkpoint
        void                                    printStoneAssignmentToWorkers(void);            //!< Prints a table showing which stones are assigned to which worker
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
        
        bool                                    resume_from_checkpoint;                         //!< Set true by initializeFromCheckpoint()
        path                                    ckp_burnin_file;                                //!< Checkpoint file path for burnin
        std::map<size_t, path>                  ckp_stone_file;                                 //!< Checkpoint file paths for stones that have one (all listed stones for flat init; first-of-sequence only for nested init)
        std::vector<std::vector<size_t>>        resume_stone_sequences;                         //!< If empty, assign stones by the usual block formula; if not, one inner vector per MPI worker
    };
    
}

#endif
