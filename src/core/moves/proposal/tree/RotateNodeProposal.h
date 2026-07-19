#ifndef RotateNodeProposal_H
#define RotateNodeProposal_H

#include <cstddef>
#include <iosfwd>
#include <vector>

#include "Proposal.h"
#include "StochasticNode.h"
#include "Tree.h"

namespace RevBayesCore {

    class TopologyNode;

    /**
     * Rotate a random internal node: permute its children, leaving every age and the clade set
     * alone.
     *
     * Child order carries no meaning in an unlabelled tree, so under most distributions this
     * proposes the state it is already in and always accepts. It is meant for a process that
     * reads the order as state: dnFBDSP takes a node's first child as the lineage that continues
     * its ancestor's species and the rest as budding descendants, so rotating a node is a
     * different budding history over the same topology. Topology moves reach those states only
     * as a side effect of changing the tree.
     *
     * The permutation is drawn uniformly over the orders that differ from the current one, in
     * both directions, so the proposal is symmetric. On a bifurcating node that is the swap.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (June Walker)
     * @since 2026-07-19, version 1.0
     *
     */
    class RotateNodeProposal : public Proposal {

    public:
        RotateNodeProposal( StochasticNode<Tree> *n );                                                //!< constructor

        // Basic utility functions
        void                                    cleanProposal(void);                                            //!< Clean up proposal
        RotateNodeProposal*                     clone(void) const;                                              //!< Clone object
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

        void                                    setChildren(size_t node_i, const std::vector<TopologyNode*> &c); //!< Rebuild a node's child list in the given order

        StochasticNode<Tree>*                   tree;                                                           //!< The variable the Proposal is working on

        size_t                                  node_index;
        std::vector<TopologyNode*>              stored_children;
        bool                                    failed;
    };

}

#endif
