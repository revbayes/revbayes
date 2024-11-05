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
        void                                                addFileMonitorExtension(const std::string &s, bool dir);
        void                                                addMonitor(const Monitor &m);
        void                                                disableScreenMonitor(bool all, size_t rep);                                             //!< Disable/remove all screen monitors
        Mcmc*                                               clone(void) const;
        void                                                checkpoint(void) const;
        void                                                finishMonitors(size_t n, MonteCarloAnalysisOptions::TraceCombinationTypes ct);          //!< Finish the monitors
        double                                              getChainLikelihoodHeat(void) const;                                                     //!< Get the heat for this chain
        double                                              getChainPosteriorHeat(void) const;                                                      //!< Get the heat for this chain
        double                                              getChainPriorHeat(void) const;
        size_t                                              getChainIndex(void) const;                                                              //!< Get the index of this chain
        const Model&                                        getModel(void) const;
        double                                              getModelLnProbability(bool like_only);
        RbVector<Monitor>&                                  getMonitors(void);
        RbVector<Move>&                                     getMoves(void);
        std::vector<tuningInfo>                             getMovesTuningInfo(void);
        MoveSchedule&                                       getSchedule(void);
        const MoveSchedule&                                 getSchedule(void) const;
        const std::string&                                  getScheduleType(void) const;
        std::string                                         getStrategyDescription(void) const;                                                     //!< Get the description of the strategy used here.
        void                                                initializeSampler(bool priorOnly=false);                                                //!< Initialize objects for mcmc sampling
        void                                                initializeSamplerFromCheckpoint( void );                                                //!< Initialize the MCMC sampler form the checkpoint file.
        void                                                monitor(std::uint64_t g);
        void                                                nextCycle(bool advanceCycle);
        bool                                                isChainActive(void);
        void                                                printOperatorSummary(bool current_period);
        void                                                redrawStartingValues(void);                                                             //!< Redraw the starting values.
        void                                                removeMonitors(void);
        void                                                reset(void);                                                                            //!< Reset the sampler and set all the counters back to 0.
        void                                                setChainActive(bool tf);
        void                                                setChainLikelihoodHeat(double v);                                                       //!< Set the heating temparature of the likelihood of the chain
        void                                                setChainPosteriorHeat(double v);                                                        //!< Set the heating temparature of the posterior of the chain
        void                                                setChainPriorHeat(double v);
        void                                                setChainIndex(size_t idx);                                                              //!< Set the index of the chain
        void                                                setCheckpointFile(const path &f);
        void                                                setLikelihoodHeat(double v);                                                            //!< Set the heating temparature of the likelihood of the chain
        void                                                setModel(Model *m, bool redraw);
        void                                                setMoves(const RbVector<Move> &mvs);
        void                                                setMovesTuningInfo(const std::vector<tuningInfo> &mvs_ti);
        void                                                setScheduleType(const std::string &s);                                                  //!< Set the type of the move schedule
        void                                                startMonitors(size_t numCycles, bool reopen);                                           //!< Start the monitors
        void                                                tune(void);                                                                             //!< Tune the sampler and its moves.
        void                                                writeMonitorHeaders(bool screen_only);                                                              //!< Write the headers of the monitors
        
        
    protected:
        void                                                resetVariableDagNodes(void);                                                //!< Extract the variable to be monitored again.
        void                                                initializeMonitors(void);                                                               //!< Assign model and mcmc ptrs to monitors
        void                                                replaceDag(const RbVector<Move> &mvs, const RbVector<Monitor> &mons);
        void                                                setActivePIDSpecialized(size_t a, size_t n);                                            //!< Set the number of processes for this class.

        
        bool                                                chain_active;
        double                                              chain_likelihood_heat;
        double                                              chain_posterior_heat;
        double                                              chain_prior_heat;
        size_t                                              chain_idx;
        path                                                checkpoint_file_name;
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
