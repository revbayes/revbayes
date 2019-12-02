#ifndef LevyJumpProposal_H
#define LevyJumpProposal_H

#include <set>
#include <string>

#include "Proposal.h"
#include "ContinuousStochasticNode.h"

namespace RevBayesCore {
    
    /**
     * The sliding operator.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Michael Landis)
     * @since 2009-09-08, version 1.0
     *
     */
    class LevyJumpProposal : public Proposal {
        
    public:
        LevyJumpProposal( ContinuousStochasticNode *n, double d);                                                                   //!<  constructor
        
        // Basic utility functions
        void                                    cleanProposal(void);                                                                //!< Clean up proposal
        LevyJumpProposal*                       clone(void) const;                                                                  //!< Clone object
        double                                  doProposal(void);                                                                   //!< Perform proposal
        const std::string&                      getProposalName(void) const;                                                        //!< Get the name of the proposal for summary printing
        double                                  getProposalTuningParameter(void) const;
        void                                    printParameterSummary(std::ostream &o, bool name_only) const;                       //!< Print the parameter summary
        void                                    prepareProposal(void);                                                              //!< Prepare the proposal
        void                                    setProposalTuningParameter(double tp);
        void                                    tune(double r);                                                                     //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                    undoProposal(void);                                                                 //!< Reject the proposal
        
    protected:
        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);                                     //!< Swap the DAG nodes on which the Proposal is working on
        
    private:
        // parameters
        
        ContinuousStochasticNode*               variable;                                                                           //!< The variable the Proposal is working on
        double                                  storedValue;                                                                        //!< The stored value of the Proposal used for rejections.
        double                                  delta;                                                                              //!< The value we propose.
    };
    
}

#endif

