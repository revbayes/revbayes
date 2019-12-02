#ifndef UnconstrainedSBNMetaMove_H
#define UnconstrainedSBNMetaMove_H

#include "AbstractMove.h"
#include "UnrootedUnconstrainedSBNProposal.h"
#include "NearestNeighborInterchange_nonClockProposal.h"
#include "SubtreePruneRegraftProposal.h"
#include "Tree.h"

#include <set>
#include <vector>

namespace RevBayesCore {

    class Proposal;

    /**
     * This move couples a SBN move to an NNI and an SPR, allowing the SBN to obtain tree samples any time a tree move is called.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-03-26, version 1.0
     *
     */
    class UnconstrainedSBNMetaMove : public AbstractMove {

    public:
        UnconstrainedSBNMetaMove(UnrootedUnconstrainedSBNProposal* u, Proposal* n, Proposal* s, std::vector<double> m, size_t n0, size_t c0, double w, bool autoTune = false);                                           //!< Constructor
        UnconstrainedSBNMetaMove(const UnconstrainedSBNMetaMove &m);                                                        //!< Copy constructor
        virtual                                                ~UnconstrainedSBNMetaMove(void);                           //!< Destructor

        // overloaded operators
        UnconstrainedSBNMetaMove&                                 operator=(const UnconstrainedSBNMetaMove &m);             //!< Assignment operator

        // pure virtual public methods
        virtual UnconstrainedSBNMetaMove*                         clone(void) const;
        size_t                                                  getActiveProposalIndex(void) const;
        const std::string&                                      getMoveName(void) const;                                //!< Get the name of the move for summary printing
        double                                                  getMoveTuningParameter(void) const;
        size_t                                                  getNumberAcceptedCurrentPeriod(void) const;             //!< Get update weight of InferenceMove
        size_t                                                  getNumberAcceptedTotal(void) const;                     //!< Get update weight of InferenceMove
        // Proposal&                                               getProposal(void);                                      //!< Get the proposal of the move
        void                                                    printSummary(std::ostream &o, bool current_period) const;                    //!< Print the move summary
        void                                                    recordState(void);
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
        size_t                                                  num_accepted_current_period;                            //!< Number of times accepted
        size_t                                                  num_accepted_total;                                     //!< Number of times accepted
        Proposal*                                               sbn;
        Proposal*                                               nni;
        Proposal*                                               spr;
        std::vector<unsigned int>                               per_move_num_accepted_current_period;
        std::vector<unsigned int>                               per_move_num_tried_current_period;
        size_t                                                  waitBeforeRecording;                                     //!< The number of moves to wait before we record any samples
        size_t                                                  waitBeforeFitting;                              //!< The number of moves to wait before we record any samples
        std::vector<Tree>                                       trace;                                                   //!< The recorded samples
        std::vector<double>                                     relativeMoveWeights;                                     //!< The relative frequency with which the moves will be used
    };

}


#endif
