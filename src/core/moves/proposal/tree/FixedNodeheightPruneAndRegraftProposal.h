#ifndef FixedNodeheightPruneAndRegraftProposal_H
#define FixedNodeheightPruneAndRegraftProposal_H

#include <cstddef>
#include <iosfwd>
#include <vector>

#include "RbVector.h"
#include "Proposal.h"

namespace RevBayesCore {

    class DagNode;
    class TopologyNode;
    class Tree;
    template <class variableType> class StochasticNode;
    
    /**
     * The fixed node-height prune-and-regraft operator.
     *
     * A fixed node-height prune-and-regraft proposal is a SPR (subtree prune-and-regraft) proposal on rooted trees without changing the node age.
     * That is, we pick a random node which is not the root.
     * Then, we prune this node and try to attach it anywhere else in the tree at this node age.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2012-07-12, version 1.0
     *
     */
    class FixedNodeheightPruneAndRegraftProposal : public Proposal {
        
    public:
        FixedNodeheightPruneAndRegraftProposal(StochasticNode<Tree> *n, StochasticNode< RbVector<Tree> > *vec_n);   //!<  constructor
        
        // Basic utility functions
        void                                    cleanProposal(void);                                                //!< Clean up proposal
        FixedNodeheightPruneAndRegraftProposal* clone(void) const;                                                  //!< Clone object
        double                                  doProposal(void);                                                   //!< Perform proposal
        const std::string&                      getProposalName(void) const;                                        //!< Get the name of the proposal for summary printing
        double                                  getProposalTuningParameter(void) const;
        void                                    prepareProposal(void);                                              //!< Prepare the proposal
        void                                    printParameterSummary(std::ostream &o, bool name_only) const;       //!< Print the parameter summary
        void                                    setProposalTuningParameter(double tp);
        void                                    tune(double r);                                                     //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                    undoProposal(void);                                                 //!< Reject the proposal
        
    protected:
        
        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);                     //!< Swap the DAG nodes on which the Proposal is working on
        
        
    private:
        // private helper methods
        void                                    findNewBrothers(std::vector<TopologyNode*> &b, TopologyNode &p, TopologyNode *n);

        // parameters
        StochasticNode<Tree>*                   variable;                                                           //!< The variable the Proposal is working on
        StochasticNode< RbVector<Tree> >*       vector_variable;                                                    //!< The alternative variable the proposal is working on

        // stored objects to undo proposal
        bool                                    failed;
        size_t                                  tree_index;
        TopologyNode*                           storedBrother;
        TopologyNode*                           storedNewBrother;

        size_t                                  storedBrotherPos;
        size_t                                  storedNewBrotherPos;
        size_t                                  storedParentPos;
    };
    
}

#endif

