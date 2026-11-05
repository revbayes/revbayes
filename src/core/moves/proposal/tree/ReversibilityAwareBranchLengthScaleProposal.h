#ifndef ReversibilityAwareBranchLengthScaleProposal_H
#define ReversibilityAwareBranchLengthScaleProposal_H

#include <cstddef>
#include <iosfwd>

#include "Proposal.h"

namespace RevBayesCore {
class DagNode;
class RateGenerator;
class Tree;
template <class variableType> class TypedDagNode;
template <class variableType> class StochasticNode;

    /**
     * A branch-length scaling proposal that treats the two branches descending from
     * the root as a single edge.
     *
     * A branch is picked at random and scaled. If that branch descends from the
     * root, its sibling is scaled by the same factor, so the pair moves together
     * and the fraction of the root edge lying on either side is left alone.
     *
     * The point of the pairing is the time-reversible model, in which the root
     * position is not identifiable and the likelihood sees the two root branches
     * only through their sum. Scaling them together moves along that identified
     * direction rather than across it.
     *
     * The pairing is nonetheless applied unconditionally, and NOT only when the
     * rate matrix is currently time reversible. Two reasons. It is a valid
     * proposal either way: scaling two branches by a common factor has Jacobian
     * sf^2, which is what the Hastings ratio below uses, and nothing about that
     * argument mentions reversibility. And making it conditional required holding
     * the rate matrix, which was expensive and fragile: a move must register every
     * node it holds, or AbstractMove::swapNode will never reach it and the pointer
     * will still refer to the original model once the model is cloned for a
     * Metropolis-coupled chain or an extra run; but registering a node also makes
     * MetropolisHastingsMove touch it, and touching the rate matrix dirties the
     * whole CTMC, turning every branch-length proposal into a full-tree likelihood
     * recomputation. Dropping the dependency avoids both.
     *
     * Note that this move must be paired with an ordinary single-branch scaler,
     * such as mvBranchLengthScale. On its own it never changes the ratio of the
     * two root branches, which under the time-reversible model is identified by
     * the prior alone and still has to be sampled.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2012-07-12, version 1.0
     *
     */
    class ReversibilityAwareBranchLengthScaleProposal : public Proposal {
        
    public:
        ReversibilityAwareBranchLengthScaleProposal( StochasticNode<Tree> *t, double d );                                  //!<  constructor
        
        // Basic utility functions
        void                                            cleanProposal(void);                                        //!< Clean up proposal
        ReversibilityAwareBranchLengthScaleProposal*    clone(void) const;                                          //!< Clone object
        double                                          doProposal(void);                                           //!< Perform proposal
        const std::string&                              getProposalName(void) const;                                //!< Get the name of the proposal for summary printing
        double                                          getProposalTuningParameter(void) const;
        void                                            prepareProposal(void);                                      //!< Prepare the proposal
        void                                            printParameterSummary(std::ostream &o, bool name_only) const;               //!< Print the parameter summary
        void                                            setProposalTuningParameter(double tp);
        void                                            tune(double r);                                             //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                            undoProposal(void);                                         //!< Reject the proposal
        
    protected:
        
        void                                            swapNodeInternal(DagNode *oldN, DagNode *newN);             //!< Swap the DAG nodes on which the Proposal is working on
        
        
    private:
        
        
        // member variables
        StochasticNode<Tree>*                           tree;

        // parameters
        double                                          delta;
                
        // stored objects to undo proposal
        double                                          stored_value;
        double                                          stored_sibling_value;
        size_t                                          stored_branch_index;
        bool                                            stored_paired;              //!< whether the sibling was scaled too, so undoProposal does not have to work it out again
    };
    
}

#endif


