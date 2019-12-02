#ifndef ElementSwapSimplexProposal_H
#define ElementSwapSimplexProposal_H

#include <set>
#include <string>

#include "SimpleProposal.h"
#include "Simplex.h"
#include "StochasticNode.h"

namespace RevBayesCore {
    
    /**
     * The element swap simplex operator.
     *
     * An element swap simplex proposal simply selects two random elements of a simplex
     * and swaps them.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2009-09-08, version 1.0
     *
     */
    class ElementSwapSimplexProposal : public SimpleProposal<Simplex> {
        
    public:
        ElementSwapSimplexProposal( StochasticNode<Simplex> *n );                                                //!<  constructor
        
        // Basic utility functions
        void                                    cleanProposal(void);                                                                //!< Clean up proposal
        ElementSwapSimplexProposal*             clone(void) const;                                                                  //!< Clone object
        double                                  propose(Simplex &v);                                                       //!< Perform proposal
        const std::string&                      getProposalName(void) const;                                                        //!< Get the name of the proposal for summary printing
        double                                  getProposalTuningParameter(void) const;
        void                                    prepareProposal(void);                                                              //!< Prepare the proposal
        void                                    printParameterSummary(std::ostream &o, bool name_only) const;                                       //!< Print the parameter summary
        void                                    setProposalTuningParameter(double tp);
        void                                    tune(double r);                                                                     //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                    undoProposal(void);                                                                 //!< Reject the proposal
        
    protected:
        
        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);                                     //!< Swap the DAG nodes on which the Proposal is working on
        
        
    private:
        // parameters
        Simplex                                 storedValue;                                                                        //!< The stored value of the Proposal used for rejections.
        //        double                              proposedValue;                                                                      //!< The value we propose.
    };
    
}

#endif
