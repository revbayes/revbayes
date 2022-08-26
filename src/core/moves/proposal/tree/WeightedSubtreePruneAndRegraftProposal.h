#ifndef WeightedSubtreePruneAndRegraftProposal_H
#define WeightedSubtreePruneAndRegraftProposal_H

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
    class WeightedSubtreePruneAndRegraftProposal : public Proposal {
        
    public:
        WeightedSubtreePruneAndRegraftProposal( StochasticNode<Tree> *n, size_t nb, double a);                      //!< constructor
        
        // Basic utility functions
        void                                    cleanProposal(void);                                                //!< Clean up proposal
        WeightedSubtreePruneAndRegraftProposal* clone(void) const;                                                  //!< Clone object
        double                                  doProposal(void);                                                   //!< Perform proposal
        const std::string&                      getProposalName(void) const;                                        //!< Get the name of the proposal for summary printing
        double                                  getProposalTuningParameter(void) const;
        void                                    prepareProposal(void);                                              //!< Prepare the proposal
        void                                    printParameterSummary(std::ostream &o, bool name_only) const;       //!< Print the parameter summary
        void                                    setProposalTuningParameter(double tp);
        void                                    tune(double r);                                                     //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                    undoProposal(void);                                                 //!< Reject the proposal
        
    protected:
        
        double                                  computeMarginal(TopologyNode& n, size_t i);                         //!< Compute the marginal over the attachement position
        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);                     //!< Swap the DAG nodes on which the Proposal is working on
        std::vector<TopologyNode*>              getSecondNodes(const TopologyNode& n);
        void                                    markDescendantNodes(const TopologyNode& n, std::vector<bool>& v);
        void                                    pruneAndRegraft(TopologyNode* a, TopologyNode* b);

        
    private:
        // private methods
        bool                                    isDescendant(const TopologyNode &n, const TopologyNode &p);

        // parameters
        StochasticNode<Tree>*                   tree;                                                               //!< The variable the Proposal is working on
        size_t                                  num_breaks;
        double                                  alpha;
        std::vector<double>                     interval;
        std::vector< std::vector<double> >      likelihoods;

        // stored objects to undo proposal
        bool                                    failed;
        TopologyNode*                           stored_choosen_node;
        TopologyNode*                           stored_sibling;
        TopologyNode*                           stored_new_sibling;
        double                                  stored_branch_length_parent;
        double                                  stored_branch_length_sibling;
        double                                  stored_branch_length_new_sibling;
    };
    
}

#endif

