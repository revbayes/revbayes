#ifndef SubtreePruneRegraftExtending_nonClockProposal_H
#define SubtreePruneRegraftExtending_nonClockProposal_H

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
    class SubtreePruneRegraftExtending_nonClockProposal : public Proposal {
        
    public:
        SubtreePruneRegraftExtending_nonClockProposal( StochasticNode<Tree> *n, double ep );                                               //!<  constructor
        
        // Basic utility functions
        void                                    cleanProposal(void);                                        //!< Clean up proposal
        SubtreePruneRegraftExtending_nonClockProposal*       clone(void) const;                                          //!< Clone object
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
        
        // parameters
        StochasticNode<Tree>*                   tree;                                                   //!< The variable the Proposal is working on
        double                                  extension_prob;
        
        // stored objects to undo proposal
        bool                                    failed;
        TopologyNode*                           storedChoosenNode;
        TopologyNode*                           storedBrother;
    };
    
}

#endif

