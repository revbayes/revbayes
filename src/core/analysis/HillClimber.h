#ifndef HillClimber_H
#define HillClimber_H

#include <cstdint>

#include "MaximumLikelihoodEstimation.h"

namespace RevBayesCore {
    
    class Monitor;
    
    /**
     * @brief Declaration of HillClimber class
     *
     * This file contains the declaration of the Markov chain Monte Carlo (HillClimber) algorithm class.
     * An HillClimber object manages the HillClimber analysis by setting up the chain, calling the moves, the monitors and etc.
     *
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2012-06-17
     *
     */
    class HillClimber : public MaximumLikelihoodEstimation {
        
    public:
        HillClimber(const Model &m, const RbVector<Move> &moves, const RbVector<Monitor> &mons);
        HillClimber(const HillClimber &m);
        virtual                                            ~HillClimber(void);                                                                             //!< Virtual destructor
        
        HillClimber&                                        operator=(const HillClimber &m);                                                               //!< Overloaded assignment operator
        
        // public methods
//        void                                                addFileMonitorExtension(const std::string &s, bool dir);
//        void                                                addMonitor(const Monitor &m);
        void                                                disableScreenMonitor(void);                                                             //!< Disable/remove all screen monitors
        HillClimber*                                        clone(void) const;
        void                                                finishMonitors(void);                                                                   //!< Finish the monitors
        Model&                                              getModel(void);
        const Model&                                        getModel(void) const;
        double                                              getModelLnProbability(bool likelihood_only);
        RbVector<Monitor>&                                  getMonitors(void);
        RbVector<Move>&                                     getMoves(void);
        MoveSchedule&                                       getSchedule(void);
        const MoveSchedule&                                 getSchedule(void) const;
        const std::string&                                  getScheduleType(void) const;
        std::string                                         getStrategyDescription(void) const;                                                     //!< Get the description of the strategy used here.
        bool                                                hasConverged(double m);
        void                                                initializeSampler(void);                                                                //!< Initialize objects for HillClimber sampling
        void                                                monitor(std::uint64_t g);
        void                                                nextCycle(void);
        void                                                removeMonitors(void);
        void                                                reset(void);                                                                            //!< Reset the sampler and set all the counters back to 0.
        void                                                setModel(Model *m);
        void                                                setNumberOfProcesses(size_t i);                                                         //!< Set the number of processes for this HillClimber simulation.
        void                                                setScheduleType(const std::string &s);                                                  //!< Set the type of the move schedule
        void                                                startMonitors(size_t num_cycles, bool reopen);                                          //!< Start the monitors
        void                                                tune(void);                                                                             //!< Tune the sampler and its moves.
        void                                                writeMonitorHeaders(void);                                                              //!< Write the headers of the monitors

    protected:
        void                                                initializeMonitors(void);                                                               //!< Assign model and HillClimber ptrs to monitors
        void                                                replaceDag(const RbVector<Move> &mvs, const RbVector<Monitor> &mons);
        void                                                setActivePIDSpecialized(size_t i, size_t n);                                                      //!< Set the number of processes for this class.
        
        
        Model*                                              model;
        RbVector<Monitor>                                   monitors;
        RbVector<Move>                                      moves;
        MoveSchedule*                                       schedule;
        std::string                                         scheduleType;                                                                           //!< Type of move schedule to be used
        
    };
    
}

#endif
