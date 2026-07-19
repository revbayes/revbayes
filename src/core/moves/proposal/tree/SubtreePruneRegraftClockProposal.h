#ifndef SubtreePruneRegraftClockProposal_H
#define SubtreePruneRegraftClockProposal_H

#include <cstddef>
#include <iosfwd>
#include <vector>

#include "Proposal.h"
#include "StochasticNode.h"
#include "Tree.h"

namespace RevBayesCore {

    class TopologyNode;

    /**
     * Prune a subtree and regraft it elsewhere, drawing a new age for the attachment.
     *
     * The clock counterpart of the unrooted subtree prune and regraft. Regrafting a subtree onto a
     * new branch of a time tree needs an age for the node that joins them, and this move draws one
     * uniformly over the ages the branch admits. Only that one age changes: every age inside the
     * moved subtree is left alone, unlike the clock nearest neighbour interchange, which rescales
     * the subtree it moves. That matters for a process whose ages are constrained by data, such as
     * the fossilized birth death range process, where a tip is an extinction and an internal node
     * is a speciation bounded by the occurrences.
     *
     * Pruning the subtree out of either state leaves the same tree, so the branch is chosen from
     * the same set in both directions and that term cancels. What is left is the node choice,
     * whose pool the regraft can change, and the two uniform age draws.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (June Walker)
     * @since 2026-07-19, version 1.0
     *
     */
    class SubtreePruneRegraftClockProposal : public Proposal {

    public:
        SubtreePruneRegraftClockProposal( StochasticNode<Tree> *n );                                  //!< constructor

        // Basic utility functions
        void                                    cleanProposal(void);                                            //!< Clean up proposal
        SubtreePruneRegraftClockProposal*       clone(void) const;                                              //!< Clone object
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

        void                                    markSubtree(const TopologyNode &n, std::vector<bool> &in) const; //!< Flag every node at or below n
        size_t                                  countMovable(void) const;                                       //!< Nodes this move may prune
        void                                    attachments(const TopologyNode &sub, std::vector<TopologyNode*> &out) const;  //!< Branches that admit the subtree
        double                                  window(const TopologyNode &sub, const TopologyNode &target) const;            //!< Ages the branch admits
        void                                    prune(TopologyNode *sub);                                       //!< Detach sub with its parent and close the gap
        void                                    attach(TopologyNode *sub, TopologyNode *target, double age);    //!< Split the branch above target and hang sub there

        StochasticNode<Tree>*                   tree;                                                           //!< The variable the Proposal is working on

        TopologyNode*                           stored_node;
        TopologyNode*                           stored_sibling;
        double                                  stored_age;
        size_t                                  stored_slot;                                                    //!< The slot the sibling had, so attach can hand it to the new branch
        bool                                    failed;
    };

}

#endif
