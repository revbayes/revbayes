#ifndef GibbsDrawCharacterHistoryProposal_H
#define GibbsDrawCharacterHistoryProposal_H

#include <set>
#include <string>

#include "StateDependentSpeciationExtinctionProcess.h"
#include "Proposal.h"
#include "StochasticNode.h"

namespace RevBayesCore {
    
    class GibbsDrawCharacterHistoryProposal : public Proposal {
        
    public:
        
        GibbsDrawCharacterHistoryProposal(StochasticNode<Tree> *n);                                                                              //!<  constructor
        void                                            cleanProposal(void);                                                                //!< Clean up proposal
        GibbsDrawCharacterHistoryProposal*              clone(void) const;                                                                  //!< Clone object
        double                                          doProposal(void);                                                                   //!< Perform proposal
        const std::string&                              getProposalName(void) const;                                                        //!< Get the name of the proposal for summary printing
        double                                          getProposalTuningParameter(void) const;

        void                                            prepareProposal(void);                                                              //!< Prepare the proposal
        void                                            printParameterSummary(std::ostream &o, bool name_only) const;                                       //!< Print the parameter summary
        void                                            setProposalTuningParameter(double tp);
        void                                            tune(double r);                                                                     //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                            undoProposal(void);                                                                 //!< Reject the proposal
        
    protected:

        void                                            swapNodeInternal(DagNode *oldN, DagNode *newN);                                     //!< Swap the DAG nodes on which the Proposal is working on
        
        StateDependentSpeciationExtinctionProcess*      distribution;

    private:
        
        StochasticNode<Tree>*                           variable;                                                                           //!< The variable the Proposal is working on
        
    };
    
}

#endif

