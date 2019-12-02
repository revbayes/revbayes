#ifndef SlideBactrianProposal_H
#define SlideBactrianProposal_H

#include <set>
#include <string>

#include "Proposal.h"
#include "ContinuousStochasticNode.h"

namespace RevBayesCore {
    
    /**
     * The sliding bactrian operator.
     *
     * A sliding bactrian proposal draws a random number from a bimodal normal distribution
     * and slides the current vale by that value.
     * Sigma is the tuning parameter of the proposal to influence the size of the proposals, the SD of the normals.
     * See Yang and Rodriguez (2013) for more details.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Andrew Magee)
     * @since 2009-09-08, version 1.0
     *
     */
    class SlideBactrianProposal : public Proposal {
        
    public:
        SlideBactrianProposal( StochasticNode<double> *n, double l, double p=0.44);                                                                      //!<  constructor
        
        // Basic utility functions
        void                                    cleanProposal(void);                                                                //!< Clean up proposal
        SlideBactrianProposal*                  clone(void) const;                                                                  //!< Clone object
        double                                  doProposal(void);                                                                   //!< Perform proposal
        const std::string&                      getProposalName(void) const;                                                        //!< Get the name of the proposal for summary printing
        double                                  getProposalTuningParameter(void) const;
        void                                    printParameterSummary(std::ostream &o, bool name_only) const;                                       //!< Print the parameter summary
        void                                    prepareProposal(void);                                                              //!< Prepare the proposal
        void                                    setProposalTuningParameter(double tp);
        void                                    tune(double r);                                                                     //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                    undoProposal(void);                                                                 //!< Reject the proposal
        
    protected:
        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);                                     //!< Swap the DAG nodes on which the Proposal is working on
        
    private:
        // parameters
        
        StochasticNode<double>*                 variable;                                                                           //!< The variable the Proposal is working on
        double                                  storedValue;                                                                        //!< The stored value of the Proposal used for rejections.
        double                                  lambda;                                                      //!< The value we propose.
    };
    
}

#endif

