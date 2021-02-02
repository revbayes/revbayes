#ifndef GammaScaleProposal_H
#define GammaScaleProposal_H

#include <iosfwd>

#include "SimpleProposal.h"

namespace RevBayesCore {
class DagNode;
template <class variableType> class StochasticNode;
    
    /**
     * The gamma scaling operator.
     *
     * A gamma scaling proposal draws a random number from a gamma distribution u ~ Gamma(lambda*1,1)
     * and scales the current vale by u
     * where lambda is the tuning parameter of the Proposal to influence the size of the proposals.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Jeremy Brown)
     * @since 2009-09-08, version 1.0
     *
     */
    class GammaScaleProposal : public SimpleProposal<double> {
        
    public:
        GammaScaleProposal( StochasticNode<double> *n, double l, double p=0.44);                                                //!<  constructor
        
        // Basic utility functions
        void                                cleanProposal(void);                                                                //!< Clean up proposal
        GammaScaleProposal*                 clone(void) const;                                                                  //!< Clone object
        const std::string&                  getProposalName(void) const;                                                        //!< Get the name of the proposal for summary printing
        double                              getProposalTuningParameter(void) const;
        void                                prepareProposal(void);                                                              //!< Prepare the proposal
        void                                printParameterSummary(std::ostream &o, bool name_only) const;                                       //!< Print the parameter summary
        double                              propose(double &v);                                                                 //!< Perform proposal
        void                                setProposalTuningParameter(double tp);
        void                                tune(double r);                                                                     //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                undoProposal(void);                                                                 //!< Reject the proposal
        
    protected:
        
        void                                swapNodeInternal(DagNode *oldN, DagNode *newN);                                     //!< Swap the DAG nodes on which the Proposal is working on
        
        
    private:
        // parameters
        
        double                              storedValue;                                                                        //!< The stored value of the Proposal used for rejections.
        double                              lambda;                                                                             //!< The scaling parameter of the Proposal
        
    };
    
}

#endif

