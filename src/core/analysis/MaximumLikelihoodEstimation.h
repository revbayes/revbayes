#ifndef MaximumLikelihoodEstimation_H
#define MaximumLikelihoodEstimation_H

#include "Cloneable.h"
#include "Model.h"
#include "Monitor.h"
#include "Move.h"
#include "MoveSchedule.h"
#include "Parallelizable.h"
#include "RandomMoveSchedule.h"
#include "RbVector.h"
#include <vector>
#include <cstdint>
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
    class MaximumLikelihoodEstimation : public Cloneable, public Parallelizable {
        
    public:
        MaximumLikelihoodEstimation(void);
        virtual                                ~MaximumLikelihoodEstimation(void);                            //!< Virtual destructor
        
        // pure virtual public methods
//        virtual void                            addFileMonitorExtension(const std::string &s, bool dir) = 0;
//        virtual void                            addMonitor(const Monitor &m) = 0;
        virtual void                            disableScreenMonitor(void) = 0;                     //!< Disable/remove all screen monitors
        virtual MaximumLikelihoodEstimation*    clone(void) const = 0;
        virtual void                            finishMonitors(void) = 0;                           //!< Finish the monitors
        virtual Model&                          getModel(void) = 0;
        virtual const Model&                    getModel(void) const = 0;
        virtual double                          getModelLnProbability(bool likelihood_only) = 0;
        virtual std::string                     getStrategyDescription(void) const = 0;             //!< Get the discription of the strategy used for this sampler.
        virtual bool                            hasConverged(double m) = 0;                         //!< Has the estimator converged to the maximum likelihood value
        virtual void                            initializeSampler(void) = 0;                        //!< Initialize objects for mcmc sampling
        virtual void                            monitor(std::uint64_t g) = 0;
        virtual void                            nextCycle(void) = 0;
        virtual void                            removeMonitors(void) = 0;
        virtual void                            reset(void) = 0;                                    //!< Reset the sampler for a new run.
        virtual void                            setModel(Model *m) = 0;
        virtual void                            startMonitors(size_t num_cycles, bool reopen) = 0;  //!< Start the monitors
        virtual void                            tune(void) = 0;                                     //!< Tune the sampler and its moves.
        virtual void                            writeMonitorHeaders(void) = 0;                      //!< Write the headers of the monitors

        // public methods
        size_t                                  getCurrentGeneration(void) const;                   //!< Get the current generations number
        //        void                                    initializeMonitors(void);                         //!< Assign model and mcmc ptrs to monitors
        //        void                                    redrawChainState(void);
        
    protected:
        
        // members
        std::uint64_t                           generation;
        
    };
    
    // Global functions using the class
    std::ostream&                               operator<<(std::ostream& o, const MaximumLikelihoodEstimation& x);                                //!< Overloaded output operator
    
}

#endif
