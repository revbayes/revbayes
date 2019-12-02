#ifndef MetropolisHastingsMetaMove_H
#define MetropolisHastingsMetaMove_H

#include "AbstractMove.h"

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
    class MetropolisHastingsMetaMove : public AbstractMove {

    public:
        MetropolisHastingsMetaMove(std::vector<Proposal*> *p, double w, bool autoTune = false);                                           //!< Constructor
        MetropolisHastingsMetaMove(const MetropolisHastingsMetaMove &m);                                                        //!< Copy constructor
        virtual                                                ~MetropolisHastingsMetaMove(void);                           //!< Destructor

        // overloaded operators
        MetropolisHastingsMetaMove&                                 operator=(const MetropolisHastingsMetaMove &m);             //!< Assignment operator

        // pure virtual public methods
        virtual MetropolisHastingsMetaMove*                         clone(void) const;
        size_t                                                  getActiveProposalIndex(void) const;
        const std::string&                                      getMoveName(void) const;                                //!< Get the name of the move for summary printing
        double                                                  getMoveTuningParameter(void) const;
        size_t                                                  getNumberAcceptedCurrentPeriod(void) const;             //!< Get update weight of InferenceMove
        size_t                                                  getNumberAcceptedTotal(void) const;                     //!< Get update weight of InferenceMove
        // Proposal&                                               getProposal(void);                                      //!< Get the proposal of the move
        void                                                    printSummary(std::ostream &o, bool current_period) const;                    //!< Print the move summary
        void                                                    setMoveTuningParameter(double tp);
        void                                                    setNumberAcceptedCurrentPeriod(size_t na);
        void                                                    setNumberAcceptedTotal(size_t na);
        void                                                    tune(void);                                             //!< Specific tuning of the move

    protected:

        //protected methods that are overwritten from the base class
        void                                                    performMcmcMove(double prHeat, double lHeat, double pHeat);            //!< Perform the move.
        void                                                    performHillClimbingMove(double lHeat, double pHeat);    //!< Perform the move.
        void                                                    resetMoveCounters(void);                                //!< Reset the counters such as numAccepted.
        virtual void                                            swapNodeInternal(DagNode *oldN, DagNode *newN);         //!< Swap the pointers to the variable on which the move works on.

    private:

        // parameters
        unsigned int                                            num_accepted_current_period;                            //!< Number of times accepted
        unsigned int                                            num_accepted_total;                                     //!< Number of times accepted
        std::vector<Proposal*>                                  proposals;                                               //!< The proposal distributions
        std::vector<unsigned int>                               per_move_num_accepted_current_period;
        std::vector<unsigned int>                               per_move_num_tried_current_period;

    };

}


#endif
