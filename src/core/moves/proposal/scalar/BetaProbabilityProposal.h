#ifndef BetaProbabilityProposal_H
#define BetaProbabilityProposal_H

#include <iosfwd>

#include "SimpleProposal.h"

namespace RevBayesCore {
class DagNode;
template <class variableType> class StochasticNode;
    
    /**
     * The beta-Probability operator.
     *
     * A beta-Probability proposal proposes a new value from a Beta distribution.
     * This proposal is applicable for parameters defined between 0 and 1.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2016-04-12, version 1.0
     *
     */
    class BetaProbabilityProposal : public SimpleProposal<double> {
        
    public:
        BetaProbabilityProposal( StochasticNode<double> *n, double d, double o);                                                                    //!<  constructor
        
        // Basic utility functions
        void                                    cleanProposal(void);                                                                //!< Clean up proposal
        BetaProbabilityProposal*                clone(void) const;                                                                  //!< Clone object
        double                                  propose(double &v);                                                                   //!< Perform proposal
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
        
        double                                  stored_value;                                                                        //!< The stored value of the Proposal used for rejections.
        double                                  delta;                                                                             //!< The scaling parameter of the Proposal
        double                                  offset;                                                                             //!< The scaling parameter of the Proposal

    };
    
}

#endif

