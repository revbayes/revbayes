#ifndef EllipticalSliceSamplingSimpleMove_H
#define EllipticalSliceSamplingSimpleMove_H

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
    class EllipticalSliceSamplingSimpleMove : public AbstractMove {

    public:
        EllipticalSliceSamplingSimpleMove(std::vector< StochasticNode<double> *> n, double window_, double weight_, size_t d, bool autoTune = false, bool forceAccept = false);        //!< Constructor
        virtual                                                 ~EllipticalSliceSamplingSimpleMove(void);                           //!< Destructor

        // public methods
        virtual EllipticalSliceSamplingSimpleMove*              clone(void) const;
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
        std::vector< StochasticNode<double> *>                  variables;                                          //!< The (Multivariate-Normally-distributed) variables the Proposal is working on
        double                                                  window;                                             //!< Window width for slice sampling
        double                                                  total_movement;                                     //!< total distance moved under auto-tuning
        int                                                     numPr;                                              //!< Number of probability evaluations
        bool                                                    force_accept;                                       //!< Should we force an acceptance of a step of size 0 after we max out on slice iterations?
    };
}


#endif
