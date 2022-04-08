#ifndef SubtreeSwapExtending_nonClockProposal_H
#define SubtreeSwapExtending_nonClockProposal_H

#include <iosfwd>

#include "Proposal.h"

namespace RevBayesCore {
class DagNode;
class TopologyNode;
class Tree;
template <class variableType> class StochasticNode;
    
    /**
     * The subtree-swap operator.
     *
     * A subtree-swap proposal is a Subtree Swap proposal on unrooted trees without changing the branch lengths.
     * That is, we randomly pick two subtrees (each defined by its MRCA node; the two nodes can not be sister, nor can they be descendant of each other) to swap.
     * When one node is an uncle of the other, then the subtree-swap proposal collapses to an NNI proposal.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2012-07-12, version 1.0
     *
     */
    class SubtreeSwapExtending_nonClockProposal : public Proposal {
        
    public:
        SubtreeSwapExtending_nonClockProposal( StochasticNode<Tree> *n, double ep );                                               //!<  constructor
        
        // Basic utility functions
        void                                                cleanProposal(void);                                        //!< Clean up proposal
        SubtreeSwapExtending_nonClockProposal*              clone(void) const;                                          //!< Clone object
        double                                              doProposal(void);                                           //!< Perform proposal
        const std::string&                                  getProposalName(void) const;                                //!< Get the name of the proposal for summary printing
        double                                              getProposalTuningParameter(void) const;
        void                                                prepareProposal(void);                                      //!< Prepare the proposal
        void                                                printParameterSummary(std::ostream &o, bool name_only) const;               //!< Print the parameter summary
        void                                                setProposalTuningParameter(double tp);
        void                                                tune(double r);                                             //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                                undoProposal(void);                                         //!< Reject the proposal
        
    protected:
        
        void                                                swapNodeInternal(DagNode *oldN, DagNode *newN);             //!< Swap the DAG nodes on which the Proposal is working on
        
        
    private:
        
        // private helper methods
        void                                                findSwappableNodes(std::vector<TopologyNode*> &b, TopologyNode &p, TopologyNode *n);
        bool                                                isDescendant(const TopologyNode &n, const TopologyNode &p);
        
        // member variables
        StochasticNode<Tree>*                               tree;
        double                                              extension_prob;
        
        // stored objects to undo proposal
        bool                                                failed;
        TopologyNode*                                       storedNodeA;
        TopologyNode*                                       storedNodeB;
        
    };
    
}

#endif

