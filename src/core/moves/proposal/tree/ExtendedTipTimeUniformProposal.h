#ifndef ExtendedTipTimeUniformProposal_H
#define ExtendedTipTimeUniformProposal_H

#include <cstddef>
#include <iosfwd>

#include "Proposal.h"
#include "StochasticNode.h"
#include "Tree.h"

namespace RevBayesCore {

    /**
     * Uniform proposal on the extinction time of a random extinct tip of an extended tree.
     *
     * The tip of an extended tree is an extinction rather than an occurrence, so it is drawn on
     * the data-fixed window (present, y_i]. A window that also tracked the augmented age would
     * depend on the current state and would need a Hastings term; the density rejects instead.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (June Walker)
     * @since 2026-07-18, version 1.0
     *
     */
    class ExtendedTipTimeUniformProposal : public Proposal {

    public:
        ExtendedTipTimeUniformProposal( StochasticNode<Tree> *n );                                       //!< constructor

        // Basic utility functions
        void                                    cleanProposal(void);                                            //!< Clean up proposal
        ExtendedTipTimeUniformProposal*         clone(void) const;                                              //!< Clone object
        double                                  doProposal(void);                                               //!< Perform proposal
        const std::string&                      getProposalName(void) const;                                    //!< Get the name of the proposal for summary printing
        double                                  getProposalTuningParameter(void) const;
        void                                    prepareProposal(void);                                          //!< Prepare the proposal
        void                                    printParameterSummary(std::ostream &o, bool name_only) const;   //!< Print the parameter summary
        void                                    setProposalTuningParameter(double tp);
        void                                    tune(double r);                                                 //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                    undoProposal(void);                                             //!< Reject the proposal

    protected:

        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);                 //!< Swap the DAG nodes on which the Proposal is working on

    private:

        StochasticNode<Tree>*                   tree;                                                           //!< The variable the Proposal is working on

        size_t                                  node_index;
        double                                  stored_age;
        bool                                    failed;
    };

}

#endif
