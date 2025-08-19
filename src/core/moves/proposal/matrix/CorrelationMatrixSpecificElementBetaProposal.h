#ifndef CorrelationMatrixSpecificElementBetaProposal_H
#define CorrelationMatrixSpecificElementBetaProposal_H

#include <cstddef>
#include <iosfwd>

#include "Proposal.h"

namespace RevBayesCore {
class DagNode;
class MatrixReal;
template <class variableType> class StochasticNode;
    
    /**
     * The Beta operator.
     *
     * This move randomly picks a specific element of a matrix of positive real numbers.
     * That means, that we select the i-th row and j-th column.
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
    class CorrelationMatrixSpecificElementBetaProposal : public Proposal {
        
    public:
        CorrelationMatrixSpecificElementBetaProposal( StochasticNode<MatrixReal> *n, size_t i, size_t j, double a, double p=0.234);                                                                      //!<  constructor
        
        // Basic utility functions
        void                                     cleanProposal(void);                                                                //!< Clean up proposal
        CorrelationMatrixSpecificElementBetaProposal*    clone(void) const;                                                                  //!< Clone object
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
        
        double                                   alpha;                                                                             //!< The Beta parameter of the move (larger lambda -> larger proposals).
        //!< The two indices of the last modified SpecificElement.
        size_t                                   indexa;
        size_t                                   indexb;
        double                                   storedValue;                                                                          //!< The value we propose.
    
    };
    
}

#endif

