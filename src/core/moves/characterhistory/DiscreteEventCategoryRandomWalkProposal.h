#ifndef DiscreteEventCategoryRandomWalkProposal_H
#define DiscreteEventCategoryRandomWalkProposal_H

#include <set>
#include <string>

#include "HeterogeneousRateBirthDeath.h"
#include "Proposal.h"
#include "StochasticNode.h"

namespace RevBayesCore {
    
    /**
     * The birth-death proposal for events along a tree.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team
     * @since 2009-09-08, version 1.0
     *
     */
    class DiscreteEventCategoryRandomWalkProposal : public Proposal {
        
    public:
        DiscreteEventCategoryRandomWalkProposal( StochasticNode<Tree> *n);                                                                //!< Constructor
        
        // Basic utility functions
        bool                                            allowClamped() const override { return true; }                                    //!< Proposal doesn't change the tree, but changes parameters describing the process that generates the tree. See #600
        void                                            cleanProposal(void) override;                                                     //!< Clean up proposal
        DiscreteEventCategoryRandomWalkProposal*        clone(void) const override;                                                       //!< Clone object
        double                                          doProposal(void) override;                                                        //!< Perform proposal
        const std::string&                              getProposalName(void) const override;                                             //!< Get the name of the proposal for summary printing
        double                                          getProposalTuningParameter(void) const override;
        void                                            printParameterSummary(std::ostream &o, bool name_only) const override;            //!< Print the parameter summary
        void                                            prepareProposal(void) override;                                                   //!< Prepare the proposal
        void                                            setProposalTuningParameter(double tp) override;
        void                                            tune(double r) override;                                                          //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                            undoProposal(void) override;                                                      //!< Reject the proposal
        
    protected:
        
        void                                            swapNodeInternal(DagNode *oldN, DagNode *newN) override;                          //!< Swap the DAG nodes on which the Proposal is working on
        
    private:
        // parameters
        
        StochasticNode<Tree>*                           variable;                                                                         //!< The variable the Proposal is working on
        HeterogeneousRateBirthDeath*                    distribution;
        
        CharacterEventDiscrete*                         stored_value;                                                                     //!< The stored value of the Proposal used for rejections.
        size_t                                          stored_category;                                                                  //!< The value we propose.
        bool                                            failed;
    };
    
}

#endif

