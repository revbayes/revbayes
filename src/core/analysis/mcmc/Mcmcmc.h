#ifndef Mcmcmc_H
#define Mcmcmc_H

#include "Mcmc.h"
#include "Model.h"
#include "Monitor.h"
#include "MonteCarloSampler.h"
#include "Move.h"

#include <vector>

namespace RevBayesCore {
    
    /**
     * @brief Parallel Metropolis-Coupled Markov chain Monte Carlo (MCMCMC) algorithm class.
     *
     * This file contains the declaration of the Markov chain Monte Carlo (MCMC) algorithm class.
     * An MCMC object manages the MCMC analysis by setting up the chain, calling the moves, the monitors and etc.
     *
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Michael Landis & Sebastian Hoehna)
     * @since Version 1.0, 2013-05-20
     *
     */
    class Mcmcmc : public MonteCarloSampler {
        
    public:
        enum class boundary
        {
            hottest,
            coldest,
            intermediate
        };

        Mcmcmc(const Model& m, const RbVector<Move> &mv, const RbVector<Monitor> &mn, std::string sT="random", size_t nc=4, size_t si=100, double dt=0.1, size_t ntries=1000, bool th=true, double tht=0.23, std::string sm="neighbor", std::string smo="multiple");
        Mcmcmc(const Mcmcmc &m);
        virtual                                ~Mcmcmc(void);                                                                   //!< Virtual destructor
        
        // public methods
        void                                    addFileMonitorExtension(const std::string &s, bool dir);
        void                                    addMonitor(const Monitor &m);
        void                                    disableScreenMonitor(bool all, size_t rep);                                     //!< Disable/remove all screen monitors
        Mcmcmc*                                 clone(void) const;
        void                                    checkpoint(void) const;
        void                                    finishMonitors(size_t n, MonteCarloAnalysisOptions::TraceCombinationTypes ct);  //!< Finish the monitors
        const Model&                            getModel(void) const;
        double                                  getModelLnProbability(bool likelihood_only);
        RbVector<Monitor>&                      getMonitors( void );
        std::string                             getStrategyDescription(void) const;                                             //!< Get the discription of the strategy used for this sampler.
        void                                    initializeSampler(bool priorOnly=false, bool suppressCharacterData=false);      //!< Initialize objects for mcmc sampling
        void                                    initializeSamplerFromCheckpoint( void );                                        //!< Initialize the MCMCMC sampler form the checkpoint file.
        void                                    monitor(unsigned long g);
        void                                    nextCycle(bool advanceCycle);
        void                                    printMoveSummary(std::ostream &o, size_t chainId, size_t moveId, Move &mv, bool current_period) const;
        void                                    printOperatorSummary(bool current_period);
        void                                    printSwapSummary(std::ostream &o) const;
        void                                    printTripSummary(std::ostream &o) const;
        void                                    printHeatSummary(std::ostream &o) const;
        void                                    printSwapSummaryPair(std::ostream &o, const size_t &row, const size_t &col) const;
        void                                    redrawStartingValues(void);                                                     //!< Redraw the starting values.
        void                                    removeMonitors(void);
        void                                    reset(void);                                                                    //!< Reset the sampler for a new run.
        void                                    resetCounters(void);                                                            //!< Reset the counters.
        void                                    setCheckpointFile(const path &f);
        void                                    setHeatsInitial(const std::vector<double> &ht);
        void                                    setSwapInterval2(const size_t &si2);
        void                                    setLikelihoodHeat(double h);                                                    //!< Set the heat of the likelihood function.
        void                                    setModel(Model *m, bool redraw);
        void                                    startMonitors(size_t numCycles, bool reopen);                                   //!< Start the monitors
        void                                    tune(void);                                                                     //!< Tune the sampler and its moves.
        void                                    writeMonitorHeaders(bool screen_only);                                                      //!< Write the headers of the monitors.

        
    protected:
        void                                    setActivePIDSpecialized(size_t i, size_t n);                                    //!< Set the number of processes for this class.

        
    private:
        void                                    initializeChains(void);
        void                                    swapChains(const std::string swap_method);
        void                                    swapMovesTuningInfo(RbVector<Move> &mvsj, RbVector<Move> &mvsk);
        void                                    swapNeighborChains(void);
        void                                    swapRandomChains(void);
        void                                    swapGivenChains(int j, int k, double lnProposalRatio = 0.0);
        void                                    updateTrips(int j);
        void                                    synchronizeValues(bool likelihood_only);
        void                                    synchronizeHeats(void);
        void                                    synchronizeTuningInfo(void);
        void                                    updateChainState(size_t j);
        double                                  computeBeta(double d, size_t i);                                                // incremental temperature schedule
        double                                  heatForChain(int i) const;
        bool                                    isColdChain(int i) const;
        int                                     heatIndexForChain(int i) const;
        int                                     chainForHeatIndex(int i) const;
        double                                  heatForIndex(int i) const;

        size_t                                  num_chains;
        std::vector<size_t>                     heat_ranks;
        std::vector<size_t>                     pid_per_chain;
        std::vector<Mcmc*>                      chains;
        std::vector<double>                     chain_values;
        std::vector<double>                     chain_heats;

        std::vector<boundary>                   chain_prev_boundary;                                // has the chain most recently visited the hottest or coldest temperature
        std::vector<int>                        chain_half_trips;                                   // how many trips has the chain made from hottest -> coldest or coldest to hottest
        std::vector<std::pair<int,int>>         heat_visitors;                                      // how many times is the chain at this temperator most recently hottest (first) or coldest (second)

        std::string                             schedule_type;
        size_t                                  current_generation;
        size_t                                  burnin_generation;
        size_t                                  swap_interval;                                      // the interval at which swaps (between neighbor chains if the swap method is neighbor or both, or between chains chosen randomly if the swap method is random) will be attempted.
        size_t                                  swap_interval2;                                     // the interval at which swaps between randomly chosen chains will be attempted (if the swapMethod is specified as both; otherwise it would be the same as swapInterval if not provided).
        
        size_t                                  active_chain_index;                                 // index of coldest chain, i.e. which one samples the posterior
        double                                  delta;                                              // delta-T, temperature increment for computeBeta
        std::vector<double>                     heat_temps;                                         // chain temperature vector, could be specified directly by user at the beginning, determined by delta, or changed during the burn-in by tuning
        bool                                    tune_heat;                                          // if the temperatures would be tuned
        double                                  tune_heat_target;                                   // the tuning target on the acceptance ratio of neighbor chain swaps
        bool                                    useNeighborSwapping;
        bool                                    useRandomSwapping;
        std::string                             swap_mode;                                          // whether making a single attempt per swap interval or attempt multiple (= nchains or nchains^2 for neighbor or random swaps, respectively) times.
        
        Mcmc*                                   base_chain;
        
        unsigned long                           generation;
        std::vector< std::vector<unsigned long> > num_attempted_swaps;
        std::vector< std::vector<unsigned long> > num_accepted_swaps;
        
        std::vector< std::vector<Mcmc::tuningInfo> >  chain_moves_tuningInfo;
    };
    
}


#endif /* defined(__rb_mlandis__Mcmcmc__) */
