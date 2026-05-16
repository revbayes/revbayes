#ifndef Mcmc_H
#define Mcmc_H

#include "MonteCarloSampler.h"

namespace RevBayesCore {
    
    class Monitor;
    
    /**
     * @brief Declaration of MCMC class
     * 
     * This file contains the declaration of the Markov chain Monte Carlo (MCMC) algorithm class. 
     * An MCMC object manages the MCMC analysis by setting up the chain, calling the moves, the monitors and etc.
     *
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2012-06-17
     *
     */
    class Mcmc : public MonteCarloSampler {
    
    public:
        Mcmc(const Model &m, const RbVector<Move> &moves, const RbVector<Monitor> &mons, size_t ntries=1000);
        Mcmc(const Mcmc &m);
        virtual                                            ~Mcmc(void);                                                                             //!< Virtual destructor
        
        Mcmc&                                               operator=(const Mcmc &m);                                                               //!< Overloaded assignment operator
        
        struct tuningInfo {
            size_t                                          num_tried_current_period;
            size_t                                          num_tried_total;
            size_t                                          num_accepted_current_period;
            size_t                                          num_accepted_total;
            double                                          tuning_parameter;
        };
        
        // public methods
        void                                                addFileMonitorExtension(const std::string &s, bool dir) override;
        void                                                addMonitor(const Monitor &m) override;
        void                                                disableScreenMonitor(bool all, size_t rep) override;                                    //!< Disable/remove all screen monitors
        Mcmc*                                               clone(void) const override;
        void                                                finishMonitors(size_t n, MonteCarloAnalysisOptions::TraceCombinationTypes ct) override; //!< Finish the monitors
        double                                              getChainLikelihoodHeat(void) const;                                                     //!< Get the likelihood heat for this chain
        double                                              getChainPosteriorHeat(void) const;                                                      //!< Get the posterior heat for this chain
        double                                              getChainPriorHeat(void) const;
        size_t                                              getChainIndex(void) const;                                                              //!< Get the index of this chain
        path                                                getCheckpointFile(void) const;
        const Model&                                        getModel(void) const override;
        double                                              getModelLnProbability(bool like_only) override;
        RbVector<Monitor>&                                  getMonitors(void) override;
        RbVector<Move>&                                     getMoves(void) override;
        std::vector<tuningInfo>                             getMovesTuningInfo(void);
        MoveSchedule&                                       getSchedule(void);
        const MoveSchedule&                                 getSchedule(void) const;
        const std::string&                                  getScheduleType(void) const;
        std::string                                         getStrategyDescription(void) const override;                                            //!< Get the description of the strategy used here.
        void                                                initializeSampler() override;                                                           //!< Initialize objects for mcmc sampling
        void                                                monitor(std::uint64_t g) override;
        void                                                nextCycle(bool advanceCycle) override;
        bool                                                isChainActive(void);
        void                                                printOperatorSummary(bool current_period) override;
        void                                                redrawStartingValues(void) override;                                                    //!< Redraw the starting values.
        void                                                removeMonitors(void) override;
        void                                                reset(void) override;                                                                   //!< Reset the sampler and set all the counters back to 0.
        void                                                setChainActive(bool tf);
        void                                                setChainLikelihoodHeat(double v);                                                       //!< Set the heating temparature of the likelihood of the chain
        void                                                setChainPosteriorHeat(double v);                                                        //!< Set the heating temparature of the posterior of the chain
        void                                                setChainPriorHeat(double v);
        void                                                setChainIndex(size_t idx);                                                              //!< Set the index of the chain
        void                                                setCheckpointFile(const path &f) override;
        void                                                setLikelihoodHeat(double v) override;                                                   //!< Set the heating temparature of the likelihood
        void                                                setModel(Model *m, bool redraw) override;
        void                                                setMoves(const RbVector<Move> &mvs);
        void                                                setMovesTuningInfo(const std::vector<tuningInfo> &mvs_ti);
        void                                                setScheduleType(const std::string &s);                                                  //!< Set the type of the move schedule
        void                                                startMonitors(size_t numCycles, bool reopen) override;                                  //!< Start the monitors
        void                                                tune(void) override;                                                                    //!< Tune the sampler and its moves.
        void                                                writeMonitorHeaders(bool screen_only) override;                                         //!< Write the headers of the monitors
        
        
    protected:
        void                                                fullCheckpoint(void) override;
        void                                                fullInitializeSamplerFromCheckpoint( void ) override;                                   //!< Initialize the MCMC sampler from the checkpoint file.
        void                                                resetVariableDagNodes(void);                                                            //!< Extract the variable to be monitored again.
        void                                                initializeMonitors(void);                                                               //!< Assign model and mcmc ptrs to monitors
        void                                                replaceDag(const RbVector<Move> &mvs, const RbVector<Monitor> &mons);
        void                                                setActivePIDSpecialized(size_t a, size_t n) override;                                   //!< Set the number of processes for this class.

        
        bool                                                chain_active;
        double                                              chain_likelihood_heat;
        double                                              chain_posterior_heat;
        double                                              chain_prior_heat;
        size_t                                              chain_idx;
        Model*                                              model;
        RbVector<Monitor>                                   monitors;
        RbVector<Move>                                      moves;
        std::vector<tuningInfo>                             moves_tuningInfo;
        size_t                                              num_init_attempts;
        MoveSchedule*                                       schedule;
        std::string                                         schedule_type;                                                                           //!< Type of move schedule to be used
        std::vector<DagNode*>                               variable_nodes;

    };

}

#endif
