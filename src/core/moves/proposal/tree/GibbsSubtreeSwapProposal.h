#ifndef GibbsSubtreeSwapProposal_H
#define GibbsSubtreeSwapProposal_H

#include <iosfwd>

#include "Proposal.h"
#include "TopologyNode.h"

namespace RevBayesCore {

class DagNode;
class Tree;

template <class variableType> class StochasticNode;
    
    /**
     * The Gibbs subtree swap operator.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2022-08-24, version 1.4.0-preview
     *
     */
    class GibbsSubtreeSwapProposal : public Proposal {
        
    public:
        GibbsSubtreeSwapProposal( StochasticNode<Tree> *n);                                                     //!< Constructor
        
        // Basic utility functions
        void                                    cleanProposal(void);                                            //!< Clean up proposal
        GibbsSubtreeSwapProposal*               clone(void) const;                                              //!< Clone object
        double                                  doProposal(void);                                               //!< Perform proposal
        const std::string&                      getProposalName(void) const;                                    //!< Get the name of the proposal for summary printing
        double                                  getProposalTuningParameter(void) const;
        void                                    prepareProposal(void);                                          //!< Prepare the proposal
        void                                    printParameterSummary(std::ostream &o, bool name_only) const;   //!< Print the parameter summary
        void                                    undoProposal(void);                                             //!< Reject the proposal
        
    protected:
        
        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);                 //!< Swap the DAG nodes the Proposal is working on
        std::vector<TopologyNode*>              getSecondNodes(const TopologyNode& n);
        void                                    markAncestralNodes(const TopologyNode& n, std::vector<bool>& v);
        void                                    markDescendantNodes(const TopologyNode& n, std::vector<bool>& v);
        void                                    swapNodes(TopologyNode* a, TopologyNode* b);
        
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

