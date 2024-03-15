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
     * The narrow-exchange operator.
     *
     * A narrow-exchange proposal is a NNI (nearest neighbour interchange) proposal on rooted trees without changing the node age.
     * That is, we pick a random node which is not the root and neither its parent is the root.
     * Then, we try to exchange the picked node with it's uncle. This move will automatically fail if the uncle is older than the parent.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2012-07-12, version 1.0
     *
     */
    class ReversibilityAwareBranchLengthScaleProposal : public Proposal {
        
    public:
        ReversibilityAwareBranchLengthScaleProposal( StochasticNode<Tree> *t, TypedDagNode<RateGenerator> *q, double d );                                               //!<  constructor
        
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
        TypedDagNode<RateGenerator>*                    q;

        // parameters
        double                                          delta;
                
        // stored objects to undo proposal
        double                                          stored_value;
        double                                          stored_sibling_value;
        size_t                                          stored_branch_index;
    };
    
}

#endif


