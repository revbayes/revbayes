#ifndef RandomIntegerWalkProposal_H
#define RandomIntegerWalkProposal_H

#include <iosfwd>

#include "Proposal.h"

namespace RevBayesCore {
class DagNode;
template <class variableType> class StochasticNode;
    
    /**
     * The random-interger-walk operator.
     *
     * This is a very simple move on integer numbers that proposes with probability p = 0.5
     * to increase the current value by 1 and with probability p = 0.5 to decrease the
     * current value by 1. Thus, it is a random walk but guided by the acceptance ratio.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2009-09-08, version 1.0
     *
     */
    class RandomIntegerWalkProposal : public Proposal {
        
    public:
        RandomIntegerWalkProposal( StochasticNode<long> *n);                                                                    //!<  constructor
        
        // Basic utility functions
        void                                cleanProposal(void);                                                                //!< Clean up proposal
        RandomIntegerWalkProposal*          clone(void) const;                                                                  //!< Clone object
        double                              doProposal(void);                                                                   //!< Perform proposal
        const std::string&                  getProposalName(void) const;                                                        //!< Get the name of the proposal for summary printing
        void                                prepareProposal(void);                                                              //!< Prepare the proposal
        void                                printParameterSummary(std::ostream &o, bool name_only) const;                                       //!< Print the parameter summary
        void                                undoProposal(void);                                                                 //!< Reject the proposal
        
    protected:
        
        void                                swapNodeInternal(DagNode *oldN, DagNode *newN);                                     //!< Swap the DAG nodes on which the Proposal is working on
        
        
    private:
        // parameters
        
        StochasticNode<long>*               variable;                                                                           //!< The variable the Proposal is working on
        long                                stored_value;                                                                        //!< The stored value of the Proposal used for rejections.
    };
    
}

#endif

