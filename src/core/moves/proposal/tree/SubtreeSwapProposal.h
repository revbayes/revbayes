#ifndef SubtreeSwapProposal_H
#define SubtreeSwapProposal_H

#include <iosfwd>

#include "Proposal.h"
#include "TopologyNode.h"

namespace RevBayesCore {

class DagNode;
class Tree;

template <class variableType> class StochasticNode;
    
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
    class SubtreeSwapProposal : public Proposal {
        
    public:
        SubtreeSwapProposal( StochasticNode<Tree> *n);                                                          //!<  constructor
        
        // Basic utility functions
        void                                    cleanProposal(void);                                            //!< Clean up proposal
        SubtreeSwapProposal*                    clone(void) const;                                              //!< Clone object
        double                                  doProposal(void);                                               //!< Perform proposal
        const std::string&                      getProposalName(void) const;                                    //!< Get the name of the proposal for summary printing
        void                                    prepareProposal(void);                                          //!< Prepare the proposal
        void                                    printParameterSummary(std::ostream &o, bool name_only) const;   //!< Print the parameter summary
        void                                    undoProposal(void);                                             //!< Reject the proposal
        
    protected:
        
        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);                 //!< Swap the DAG nodes on which the Proposal is working on
        
        
    private:
        // private methods
        bool                                    isDescendant(const TopologyNode &n, const TopologyNode &p);

        // parameters
        StochasticNode<Tree>*                   tree;                                                   //!< The variable the Proposal is working on
        
        // stored objects to undo proposal
        bool                                    failed;
        TopologyNode*                           stored_first_node;
        TopologyNode*                           stored_second_node;

    };
    
}

#endif

