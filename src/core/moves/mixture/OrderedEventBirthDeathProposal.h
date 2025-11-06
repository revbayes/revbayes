#ifndef OrderedEventBirthDeathProposal_H
#define OrderedEventBirthDeathProposal_H

#include <cstddef>
#include <iosfwd>

#include "Proposal.h"

namespace RevBayesCore {
class DagNode;
class OrderedEventTimes;
template <class variableType> class StochasticNode;
    
    class OrderedEventBirthDeathProposal : public Proposal {
        
    public:
        OrderedEventBirthDeathProposal(StochasticNode<OrderedEventTimes> *n);                                                       //!<  constructor
        
        // Basic utility functions
        void                                    cleanProposal(void);                                           //!< Clean up proposal
        OrderedEventBirthDeathProposal*         clone(void) const;                                             //!< Clone object
        double                                  doProposal(void);                                              //!< Perform proposal
        const std::string&                      getProposalName(void) const;                                   //!< Get the name of the proposal for summary printing
        double                                  getProposalTuningParameter(void) const;
        void                                    prepareProposal(void);                                         //!< Prepare the proposal
        void                                    printParameterSummary(std::ostream &o, bool name_only) const;  //!< Print the parameter summary
        void                                    setProposalTuningParameter(double tp);
        void                                    tune(double r);                                                //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                    undoProposal(void);                                            //!< Reject the proposal
        
    protected:
        
        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);                //!< Swap the DAG nodes on which the Proposal is working on
        
    private:
        
        // parameters
        StochasticNode<OrderedEventTimes>*      event_var;                                                     //!< The variable the Proposal is working on
        
        // stored objects to undo proposal
        double                                  event_time;
        bool                                    was_birth;
        
    };
    
}

#endif

