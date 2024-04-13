#ifndef EventBranchTimeBetaProposal_H
#define EventBranchTimeBetaProposal_H

#include <set>
#include <string>

#include "AbstractCharacterHistoryBirthDeathProcess.h"
#include "Proposal.h"
#include "StochasticNode.h"
#include "CharacterEvent.h"

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
    class EventBranchTimeBetaProposal : public Proposal {
        
    public:
        EventBranchTimeBetaProposal( StochasticNode<Tree> *n, double d, double o);                                                                //!<  constructor
        
        // Basic utility functions
        void                                            cleanProposal(void);                                                                //!< Clean up proposal
        EventBranchTimeBetaProposal*                    clone(void) const;                                                                  //!< Clone object
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
        AbstractCharacterHistoryBirthDeathProcess*      distribution = nullptr;
        double                                          delta = -1;
        double                                          offset = -1;

        CharacterEvent*                                 stored_value = nullptr;                                                             //!< The stored value of the Proposal used for rejections.
        double                                          stored_age = -1;                                                                    //!< The value we propose.
        size_t                                          stored_branch_index = -1;                                                           //!< The value we propose.
        bool                                            failed = false;
    };
    
}

#endif

