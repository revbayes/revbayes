#ifndef UnrootedUnconstrainedSBNProposal_H
#define UnrootedUnconstrainedSBNProposal_H

#include <string>

#include "Proposal.h"
#include "StochasticNode.h"
#include "Tree.h"
#include "SBNParameters.h"
#include "UnconstrainedSBN.h"

namespace RevBayesCore {

    /**
     * The subtree-prune-and-regraft operator.
     *
     * A subtree-prune-and-regraft proposal is a SPR proposal on unrooted trees without changing the branch lengths.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2012-07-12, version 1.0
     *
     */
    class UnrootedUnconstrainedSBNProposal : public Proposal {

    public:
        UnrootedUnconstrainedSBNProposal( StochasticNode<Tree> *n );                                               //!<  constructor
        UnrootedUnconstrainedSBNProposal( StochasticNode<Tree> *n, const SBNParameters sbn );                                               //!<  constructor

        // Basic utility functions
        void                                    cleanProposal(void);                                        //!< Clean up proposal
        UnrootedUnconstrainedSBNProposal*            clone(void) const;                                          //!< Clone object
        double                                  doProposal(void);                                           //!< Perform proposal
        void                                    fitSBNToTreeSamples(std::vector<Tree> &t);
        const std::string&                      getProposalName(void) const;                                //!< Get the name of the proposal for summary printing
        double                                  getProposalTuningParameter(void) const;
        void                                    prepareProposal(void);                                      //!< Prepare the proposal
        void                                    printParameterSummary(std::ostream &o, bool name_only) const;               //!< Print the parameter summary
        void                                    setProposalTuningParameter(double tp);
        void                                    tune(double r);                                             //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                    undoProposal(void);                                         //!< Reject the proposal

    protected:

        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);             //!< Swap the DAG nodes on which the Proposal is working on


    private:
        // private methods

        // parameters
        StochasticNode<Tree>*                   tree;                                                   //!< The variable the Proposal is working on
        SBNParameters                           sbn;
        UnconstrainedSBN                        SBNDistribution;

        // stored objects to undo proposal
        Tree                                    stored_tree;
    };

}

#endif
