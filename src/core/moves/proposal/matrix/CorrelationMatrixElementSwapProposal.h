#ifndef CorrelationMatrixElementSwapProposal_H
#define CorrelationMatrixElementSwapProposal_H

#include <cstddef>
#include <iosfwd>

#include "Proposal.h"

namespace RevBayesCore {
class DagNode;
class MatrixReal;
template <class variableType> class StochasticNode;
    
    /**
     * The Swap operator.
     *
     * This move randomly picks an element of a matrix of positive real numbers.
     * That means, that we randomly pick the i-th row and j-th column with equal probability.
     * Then, we propose a Swap distance and slide the value.
     * The actual Swap distance is computed by delta = lambda * ( u - 0.5 )
     * where u ~ Uniform(0,1).
     * The proposal is thus m[i][j] += lambda * ( u - 0.5 )
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Michael R. May)
     * @since 2009-09-08, version 1.0
     *
     */
    class CorrelationMatrixElementSwapProposal : public Proposal {
        
    public:
        CorrelationMatrixElementSwapProposal( StochasticNode<MatrixReal> *n );                                                                      //!<  constructor
        
        // Basic utility functions
        void                                     cleanProposal(void);                                                                //!< Clean up proposal
        CorrelationMatrixElementSwapProposal*    clone(void) const;                                                                  //!< Clone object
        double                                   doProposal(void);                                                                   //!< Perform proposal
        const std::string&                       getProposalName(void) const;                                                        //!< Get the name of the proposal for summary printing
        double                                   getProposalTuningParameter(void) const;
        void                                     printParameterSummary(std::ostream &o, bool name_only) const;                                       //!< Print the parameter summary
        void                                     prepareProposal(void);                                                              //!< Prepare the proposal
        void                                     setProposalTuningParameter(double tp);
        void                                     tune(double r);                                                                     //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                     undoProposal(void);                                                                 //!< Reject the proposal
        
    protected:
        void                                     swapNodeInternal(DagNode *oldN, DagNode *newN);                                     //!< Swap the DAG nodes on which the Proposal is working on
        
    private:
        // parameters
        
        StochasticNode<MatrixReal >*             variable;
        
        //!< The indices of the last modified element.
        size_t                                   indexa_1;
        size_t                                   indexb_1;
        double                                   storedValue_1;                                                                          //!< The value we propose.

        size_t                                   indexa_2;
        size_t                                   indexb_2;
        double                                   storedValue_2;                                                                          //!< The value we propose.
        
    };
    
}

#endif

