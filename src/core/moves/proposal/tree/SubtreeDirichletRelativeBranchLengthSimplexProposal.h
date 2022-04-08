#ifndef SubtreeDirichletRelativeBranchLengthSimplexProposal_H
#define SubtreeDirichletRelativeBranchLengthSimplexProposal_H

#include <iosfwd>
#include <vector>

#include "Proposal.h"
#include "Simplex.h"
#include "StochasticNode.h"
#include "Tree.h"

namespace RevBayesCore {
class DagNode;
class TopologyNode;
class Tree;
template <class variableType> class StochasticNode;
    
    /**
     * The subtree-scale operator.
     *
     * A subtree-scale proposal is a scaling proposal on rooted subtrees without changing the topology.
     * That is, we pick a random node which is not the root.
     * Then, we uniformly pick an age between the parent and the oldest sampled descendant.
     * The picked subtree is then scaled to this new age.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2012-07-12, version 1.0
     *
     */
    class SubtreeDirichletRelativeBranchLengthSimplexProposal : public Proposal {
        
    public:
        SubtreeDirichletRelativeBranchLengthSimplexProposal( StochasticNode<Tree> *tr, StochasticNode<Simplex> *relbls, double a, double p=0.44);                                                     //!<  constructor
        
        // Basic utility functions
        void                                                                 cleanProposal(void);                                        //!< Clean up proposal
        SubtreeDirichletRelativeBranchLengthSimplexProposal*                 clone(void) const;                                          //!< Clone object
        double                                                               doProposal(void);                                           //!< Perform proposal
        const std::string&                                                   getProposalName(void) const;                                //!< Get the name of the proposal for summary printing
        double                                                               getProposalTuningParameter(void) const;
        void                                                                 prepareProposal(void);                                      //!< Prepare the proposal
        void                                                                 printParameterSummary(std::ostream &o, bool name_only) const;               //!< Print the parameter summary
        void                                                                 setProposalTuningParameter(double tp);
        void                                                                 tune(double r);                                             //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                                                 undoProposal(void);                                         //!< Reject the proposal
        
    protected:
        
        void                                                                 swapNodeInternal(DagNode *oldN, DagNode *newN);             //!< Swap the DAG nodes on which the Proposal is working on
        
        
    private:
        
        // parameters
        StochasticNode<Tree>*                   tree;
        StochasticNode<Simplex>*                relative_branch_lengths;
        double                                  alpha;
        
        // stored objects to undo proposal
        bool                                    failed;
        Simplex                                 stored_relative_branch_lengths;

    };
    
}

#endif

