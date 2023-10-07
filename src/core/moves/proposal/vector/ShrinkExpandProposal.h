#ifndef ShrinkExpandProposal_H
#define ShrinkExpandProposal_H

#include <cstddef>
#include <ostream>
#include <vector>

#include "Proposal.h"

namespace RevBayesCore {
class DagNode;
template <class variableType> class StochasticNode;
    
    /**
     * @brief Scaling Proposal of a all elements of a vector.
     *
     *
     * This proposal randomly scales all elements of a vector using the same scaling factor.
     * The actual scaling factor is computed by sf = exp( lambda * ( u - 0.5 ) )
     * where u ~ Uniform(0,1).
     * It generally makes more sense to apply the scaling proposal on a vector of positive
     * real numbers but technically it works on negative numbers too. However,
     * the proposal will never change the sign of the value and thus is incomplete if applied
     * to variable defined on the whole real line.
     *
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @copyright GPL version 3
     * @since 2015-05-21, version 1.0
     *
     */
    class ShrinkExpandProposal : public Proposal {
        
    public:
        ShrinkExpandProposal(std::vector<StochasticNode<double> *> n, StochasticNode<double> *s, double l);                             //!< Constructor
        
        void                                        cleanProposal(void);                                                                //!< Clean up proposal
        ShrinkExpandProposal*                       clone(void) const;                                                                  //!< Clone object
        double                                      doProposal(void);                                                                   //!< Perform proposal
        const std::string&                          getProposalName(void) const;                                                        //!< Get the name of the proposal for summary printing
        double                                      getProposalTuningParameter(void) const;
        void                                        printParameterSummary(std::ostream &o, bool name_only) const;                                       //!< Print the parameter summary
        void                                        prepareProposal(void);                                                              //!< Prepare the proposal
        void                                        setProposalTuningParameter(double tp);
        void                                        tune(double r);                                                                     //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                        undoProposal(void);                                                                 //!< Reject the proposal
        
    protected:
        
        void                                        swapNodeInternal(DagNode *oldN, DagNode *newN);                                     //!< Swap the DAG nodes on which the Proposal is working on
        
        
    private:
        // parameters
        
        std::vector<StochasticNode<double> *>       variables;
        StochasticNode<double>                     *sd;
        double                                      lambda;                                                                             //!< The scale parameter of the Proposal (larger lambda -> larger proposals).
        size_t                                      length;
        double                                      stored_scaling_factor;                                                                        //!< The stored value of the last modified element.
        double                                      stored_mean;                                                                        //!< The stored value of the last modified element.
        
        
    };
    
}

#endif

