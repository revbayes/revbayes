#ifndef MonteCarloSampler_H
#define MonteCarloSampler_H

#include "Cloneable.h"
#include "Model.h"
#include "Monitor.h"
#include "MonteCarloAnalysis.h"
#include "Move.h"
#include "MoveSchedule.h"
#include "Parallelizable.h"
#include "RandomMoveSchedule.h"
#include "RbVector.h"
#include <vector>
#include "SequentialMoveSchedule.h"

namespace RevBayesCore {
    
    class Monitor;
    
    /**
     * @brief Abstract class for all Monte Carlo samplers, such as the Markov chain Monte Carlo (MCMC) sampler and Metropolis coupled MCMC (MC^3).
     *
     * The Monte Carlo sampler base class provides the common methods for all Monte Carlo samplers, such as
     * - Markov chain Monte Carlo (MCMC)
     * - Metropolis coupled Markov chain Monte Carlo (MCMCMC or MC^3)
     * - PowerPosteriorSampler
     *
     *
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2014-06-19
     *
     */
    class MonteCarloSampler : public Cloneable, public Parallelizable {
        
    public:
        MonteCarloSampler(void);
        MonteCarloSampler(const MonteCarloSampler &m);
        virtual                                ~MonteCarloSampler(void);                                    //!< Virtual destructor
                
        // pure virtual public methods
        virtual void                            addFileMonitorExtension(const std::string &s, bool dir) = 0;
        virtual void                            addMonitor(const Monitor &m) = 0;
        virtual void                            disableScreenMonitor(bool all, size_t rep) = 0;             //!< Disable/remove all screen monitors
        virtual MonteCarloSampler*              clone(void) const = 0;
        virtual void                            checkpoint(void) const = 0;                                 //!< Perform checkpointing by writing the current values to a file.
//        virtual void                            run(size_t g) = 0;
        virtual void                            finishMonitors(size_t n, MonteCarloAnalysisOptions::TraceCombinationTypes ct) = 0; //!< Finish the monitors
        virtual const Model&                    getModel(void) const = 0;
        virtual double                          getModelLnProbability(bool like_only) = 0;
        virtual RbVector<Monitor>&              getMonitors() = 0;
        virtual std::string                     getStrategyDescription(void) const = 0;                     //!< Get the discription of the strategy used for this sampler.
        virtual void                            initializeSamplerFromCheckpoint(void) = 0;                  //!< Initialize the values from the checkpoint file
        virtual void                            initializeSampler(bool priorOnly=false) = 0;                //!< Initialize objects for mcmc sampling
        virtual void                            monitor(std::uint64_t g) = 0;
        virtual void                            nextCycle(bool advanceCycle) = 0;
        virtual void                            printOperatorSummary(bool current_period) = 0;
        virtual void                            redrawStartingValues(void) = 0;                             //!< Redraw the starting values.
        virtual void                            removeMonitors(void) = 0;
        virtual void                            reset(void) = 0;                                            //!< Reset the sampler for a new run.
        virtual void                            setCheckpointFile(const path &f) = 0;
        virtual void                            setLikelihoodHeat(double v) = 0;                            //!< Set the heating temparature of the likelihood of the chain
//        virtual void                            setMasterSampler(bool tf) = 0;                            //!< Set whether this one is the master.
        virtual void                            setModel(Model *m, bool redraw) = 0;
        virtual void                            startMonitors(size_t numCycles, bool reopen) = 0;           //!< Start the monitors
        virtual void                            tune(void) = 0;                                             //!< Tune the sampler and its moves.
        virtual void                            writeMonitorHeaders(bool screen_only) = 0;                              //!< Write the headers of the monitors

        // public methods
        size_t                                  getCurrentGeneration(void) const;                           //!< Get the current generations number
        void                                    setCurrentGeneration(size_t g);
        //        void                                    initializeMonitors(void);                         //!< Assign model and mcmc ptrs to monitors
//        void                                    redrawChainState(void);
        
    protected:
                
        // members
        std::uint64_t                           generation;
        
    };

    // Global functions using the class
    std::ostream&                               operator<<(std::ostream& o, const MonteCarloSampler& x);                                //!< Overloaded output operator

}

#endif
