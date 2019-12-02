#ifndef CorrelationMatrixElementBetaProposal_H
#define CorrelationMatrixElementBetaProposal_H

#include <set>
#include <string>

#include "Proposal.h"
#include "MatrixReal.h"
#include "StochasticNode.h"

namespace RevBayesCore {
    
    /**
     * The Beta operator.
     *
     * This move randomly picks an element of a matrix of positive real numbers.
     * That means, that we randomly pick the i-th row and j-th column with equal probability.
     * Then, we propose a Beta distance and slide the value.
     * The actual Beta distance is computed by delta = lambda * ( u - 0.5 )
     * where u ~ Uniform(0,1).
     * The proposal is thus m[i][j] += lambda * ( u - 0.5 )
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Michael R. May)
     * @since 2009-09-08, version 1.0
     *
     */
    class CorrelationMatrixRandomWalkProposal : public Proposal {
        
    public:
        CorrelationMatrixRandomWalkProposal( StochasticNode<MatrixReal> *n, double s, double p=0.44);                                                                      //!<  constructor
        
        // Basic utility functions
        void                                     cleanProposal(void);                                                                //!< Clean up proposal
        CorrelationMatrixRandomWalkProposal*     clone(void) const;                                                                  //!< Clone object
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
        double                                   sigma;
        MatrixReal                               stored_matrix;                                                                          //!< The value we propose.
        
    };
    
}

#endif

