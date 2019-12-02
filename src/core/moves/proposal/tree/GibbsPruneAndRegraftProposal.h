#ifndef GibbsPruneAndRegraftProposal_H
#define GibbsPruneAndRegraftProposal_H

#include <string>

#include "Proposal.h"
#include "StochasticNode.h"
#include "Tree.h"

namespace RevBayesCore {
    
    /**
     * The metropolized Gibbs fixed node-height prune-and-regraft operator.
     *
     * A metropolized Gibbs fixed node-height prune-and-regraft proposal is a SPR (subtree prune-and-regraft) proposal on rooted trees without changing the node age.
     * That is, we pick a random node which is not the root.
     * Then, we prune this node and try to attach it at all possible re-attachment points elsewhere in the tree at this node age.
     * Finally, we pick the re-attachment point according to the tree probability.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2012-07-12, version 1.0
     *
     */
    class GibbsPruneAndRegraftProposal : public Proposal {
        
    public:
        GibbsPruneAndRegraftProposal( StochasticNode<Tree> *n);                                             //!<  constructor
        
        // Basic utility functions
        void                                    cleanProposal(void);                                        //!< Clean up proposal
        GibbsPruneAndRegraftProposal*           clone(void) const;                                          //!< Clone object
        double                                  doProposal(void);                                           //!< Perform proposal
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
        // private helper methods
        void                                    findNewBrothers(std::vector<TopologyNode*> &b, TopologyNode &p, TopologyNode *n);
        TopologyNode*                           pruneAndRegraft(TopologyNode *brother, TopologyNode *newBrother, TopologyNode *parent, TopologyNode &grandparent);
        
        // parameters
        StochasticNode<Tree>*                   variable;                                                   //!< The variable the Proposal is working on
        
        // stored objects to undo proposal
        bool                                    failed;
        TopologyNode*                           storedBrother;
        TopologyNode*                           storedNewBrother;
    };
    
}

#endif

