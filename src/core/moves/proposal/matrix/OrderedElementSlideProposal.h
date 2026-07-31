#ifndef OrderedElementSlideProposal_H
#define OrderedElementSlideProposal_H

#include <cstddef>
#include <iosfwd>

#include "Proposal.h"

namespace RevBayesCore {
class DagNode;
class MatrixReal;
template <class variableType> class StochasticNode;

    /**
     * Slide one entry of an ordered vector, uniformly between the entries either side of it.
     *
     * The value must expose its parts as ordered vectors; the proposal asks it for one, picks a
     * free entry, and redraws that entry inside the interval its neighbours leave. Ordering is
     * therefore preserved by construction rather than repaired or rejected afterwards, and the
     * proposal is symmetric, so the Hastings ratio is one.
     *
     * There is no tuning parameter. The window is the local gap, so it is already the right size
     * wherever it is applied, in the way mvNodeTimeSlideUniform is for node ages.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (June Walker)
     * @since 2026-07-31, version 1.0
     */
    class OrderedElementSlideProposal : public Proposal {

    public:
        OrderedElementSlideProposal( StochasticNode<MatrixReal> *n );

        void                                    cleanProposal(void);
        OrderedElementSlideProposal*            clone(void) const;
        double                                  doProposal(void);
        const std::string&                      getProposalName(void) const;
        double                                  getProposalTuningParameter(void) const;
        void                                    prepareProposal(void);
        void                                    printParameterSummary(std::ostream &o, bool name_only) const;
        void                                    setProposalTuningParameter(double tp);
        void                                    tune(double r);
        void                                    undoProposal(void);

    protected:
        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);

    private:
        StochasticNode<MatrixReal>*             variable;

        size_t                                  stored_vector;
        size_t                                  stored_entry;
        double                                  stored_value;
        bool                                    failed;
    };

}

#endif
