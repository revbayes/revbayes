#ifndef EllipticalSliceSamplingLognormalIIDMove_H
#define EllipticalSliceSamplingLognormalIIDMove_H

#include "AbstractMove.h"
#include "StochasticNode.h"

#include <set>
#include <vector>

namespace RevBayesCore {
    
    class Proposal;
    
    /**
     * Base class for all Metropolis-Hastings within an MCMC. 
     *
     * The base class of all moves only provides the interface for the call to propose a new move.
     * Here the perform methods actually does the accept/reject step.
     * All specifics are implemented in the derived classes.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-03-26, version 1.0
     *
     */
    class EllipticalSliceSamplingLognormalIIDMove : public AbstractMove {

    public:
        EllipticalSliceSamplingLognormalIIDMove( std::vector< StochasticNode<double> *> n, TypedDagNode<double>* m, TypedDagNode<double>* s, double window_, double weight_, size_t d, bool t );        //!< Constructor
        virtual                                                 ~EllipticalSliceSamplingLognormalIIDMove(void);                           //!< Destructor

        // public methods
        virtual EllipticalSliceSamplingLognormalIIDMove*              clone(void) const;
        const std::string&                                      getMoveName(void) const;                            //!< Get the name of the move for summary printing
        double                                                  getMoveTuningParameter(void) const;
        void                                                    printSummary(std::ostream &o, bool current_period) const;                //!< Print the move summary
        void                                                    setMoveTuningParameter(double tp);
        void                                                    tune(void);                                         //!< Specific tuning of the move

    protected:
        //protected methods that are overwritten from the base class
        void                                                    performMcmcMove(double prHeat, double lHeat, double pHeat);            //!< Perform the move.
        void                                                    resetMoveCounters(void);                            //!< Reset the counters such as numAccepted.
        virtual void                                            swapNodeInternal(DagNode *oldN, DagNode *newN);             //!< Swap the pointers to the variable on which the move works on.
        
    private:

        // parameters
        std::vector< StochasticNode<double> *>                  variables;                                          //!< The (iid lognormal) variables the Proposal is working on
        TypedDagNode<double>*                                   mu;                                                 //!< The log-mean parameter of the lognormals
        TypedDagNode<double>*                                   sigma;                                              //!< The log-sd parameter of the lognormals
        double                                                  window;                                             //!< Window width for slice sampling
        double                                                  total_movement;                                     //!< total distance moved under auto-tuning
        int                                                     numPr;                                              //!< Number of probability evaluations
    };
}


#endif
