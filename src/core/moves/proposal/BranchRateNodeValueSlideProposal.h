#ifndef BranchRateNodeValueSlideProposal_H
#define BranchRateNodeValueSlideProposal_H

#include <iosfwd>

#include "Proposal.h"
#include "TypedDagNode.h"
#include "StochasticNode.h"
#include "RbVector.h"
#include "Tree.h"

namespace RevBayesCore {
    
    /**
     * The scaling operator.
     *
     * A scaling proposal draws a random uniform number u ~ unif (-0.5,0.5)
     * and Slides the current vale by a scaling factor
     * sf = exp( lambda * u )
     * where lambda is the tuning parameter of the Proposal to influence the size of the proposals.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2009-09-08, version 1.0
     *
     */
    class BranchRateNodeValueSlideProposal : public Proposal {
        
    public:
        BranchRateNodeValueSlideProposal( StochasticNode< RbVector<double> > *n, TypedDagNode< Tree > *t, double l, double p=0.44);                                                                    //!<  constructor
        
        // Basic utility functions
        void                                cleanProposal(void);                                                                //!< Clean up proposal
        BranchRateNodeValueSlideProposal*   clone(void) const;                                                                  //!< Clone object
        double                              doProposal(void);                                                                   //!< Perform proposal
        const std::string&                  getProposalName(void) const;                                                        //!< Get the name of the proposal for summary printing
        double                              getProposalTuningParameter(void) const;
        void                                prepareProposal(void);                                                              //!< Prepare the proposal
        void                                printParameterSummary(std::ostream &o, bool name_only) const;                                       //!< Print the parameter summary
        void                                setProposalTuningParameter(double tp);
        void                                tune(double r);                                                                     //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                undoProposal(void);                                                                 //!< Reject the proposal
        
    protected:
        
        void                                swapNodeInternal(DagNode *oldN, DagNode *newN);                                     //!< Swap the DAG nodes on which the Proposal is working on
        
        
    private:
        // parameters
        
        StochasticNode< RbVector<double> >* variable;                                                                           //!< The variable the Proposal is working on
        TypedDagNode<Tree>*                 tree;                                                                           //!< The variable the Proposal is working on
        double                              stored_value;                                                                        //!< The stored value of the Proposal used for rejections.
        size_t                              stored_index;                                                                        //!< The stored value of the Proposal used for rejections.
        double                              lambda;                                                                             //!< The scaling parameter of the Proposal

    };
    
}

#endif

