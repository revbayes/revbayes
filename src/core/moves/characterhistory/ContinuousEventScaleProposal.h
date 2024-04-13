#ifndef ContinuousEventScaleProposal_H
#define ContinuousEventScaleProposal_H

#include <set>
#include <string>

#include "ConditionedBirthDeathShiftProcessContinuous.h"
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
    class ContinuousEventScaleProposal : public Proposal {
        
    public:
        ContinuousEventScaleProposal( StochasticNode<Tree> *n, double l );                                                                //!<  constructor
        
        // Basic utility functions
        void                                            cleanProposal(void);                                                                //!< Clean up proposal
        ContinuousEventScaleProposal*                   clone(void) const;                                                                  //!< Clone object
        double                                          doProposal(void);                                                                   //!< Perform proposal
        const std::string&                              getProposalName(void) const;                                                        //!< Get the name of the proposal for summary printing
        double                                          getProposalTuningParameter(void) const;
        void                                            printParameterSummary(std::ostream &o, bool name_only) const;                                       //!< Print the parameter summary
        void                                            prepareProposal(void);                                                              //!< Prepare the proposal
        void                                            setProposalTuningParameter(double tp);
        void                                            tune(double r);                                                                     //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                            undoProposal(void);                                                                 //!< Reject the proposal
        
    protected:
        
        void                                            swapNodeInternal(DagNode *oldN, DagNode *newN);                                     //!< Swap the DAG nodes on which the Proposal is working on
        
    private:
        // parameters
        
        StochasticNode<Tree>*                           variable = nullptr;                                                                 //!< The variable the Proposal is working on
        ConditionedBirthDeathShiftProcessContinuous*    distribution = nullptr;

        double                                          lambda = -1;
        
        CharacterEventContinuous*                       stored_value = nullptr;                                                              //!< The stored value of the Proposal used for rejections.
        double                                          stored_rate = -1;                                                                    //!< The value we propose.
        size_t                                          stored_index = -1;                                                                   //!< The value we propose.
        bool                                            failed = false;

};
    
}

#endif


