#ifndef WeightedSubtreeSwapProposal_H
#define WeightedSubtreeSwapProposal_H

#include <iosfwd>

#include "Proposal.h"
#include "TopologyNode.h"

namespace RevBayesCore {

class DagNode;
class Tree;

template <class variableType> class StochasticNode;
    
    /**
     * The weighted subtree swap operator.
     *
     * A weighted subtree swap proposal jointly updates the tree topology and branch lengths.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2022-08-26, version 1.4.0-preview
     *
     */
    class WeightedSubtreeSwapProposal : public Proposal {
        
    public:
        WeightedSubtreeSwapProposal( StochasticNode<Tree> *n, size_t nb, double a);                                 //!< Constructor
        
        // Basic utility functions
        void                                    cleanProposal(void);                                                //!< Clean up proposal
        WeightedSubtreeSwapProposal*            clone(void) const;                                                  //!< Clone object
        double                                  doProposal(void);                                                   //!< Perform proposal
        const std::string&                      getProposalName(void) const;                                        //!< Get the name of the proposal for summary printing
        double                                  getProposalTuningParameter(void) const;
        void                                    prepareProposal(void);                                              //!< Prepare the proposal
        void                                    printParameterSummary(std::ostream &o, bool name_only) const;       //!< Print the parameter summary
        void                                    undoProposal(void);                                                 //!< Reject the proposal
        
    protected:
        
        double                                  computeMarginal(TopologyNode& n, size_t i);                         //!< Compute the marginal over the attachement position
        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);                     //!< Swap the DAG nodes the Proposal is working on
        std::vector<TopologyNode*>              getSecondNodes(const TopologyNode& n);
        void                                    markAncestralNodes(const TopologyNode& n, std::vector<bool>& v);
        void                                    markDescendantNodes(const TopologyNode& n, std::vector<bool>& v);
        void                                    swapNodes(TopologyNode* a, TopologyNode* b);
        
    private:
        // private methods
        bool                                    isDescendant(const TopologyNode &n, const TopologyNode &p);

        // parameters
        StochasticNode<Tree>*                   tree;                                                               //!< The variable the Proposal is working on
        size_t                                  num_breaks;
        double                                  alpha;                                                              //!< Not a tuning parameter!
        std::vector<double>                     interval;
        std::vector< std::vector<double> >      likelihoods;

        // stored objects to undo proposal
        bool                                    failed;
        TopologyNode*                           stored_first_node;
        TopologyNode*                           stored_second_node;
        double                                  stored_branch_length_parent;
        double                                  stored_branch_length_sibling;
        double                                  stored_branch_length_new_sibling;

    };
    
}

#endif

