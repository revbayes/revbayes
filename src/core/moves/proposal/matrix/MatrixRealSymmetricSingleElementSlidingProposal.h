#ifndef MatrixRealSymmetricSingleElementSlidingProposal_H
#define MatrixRealSymmetricSingleElementSlidingProposal_H

#include <cstddef>
#include <iosfwd>

#include "Proposal.h"

namespace RevBayesCore {
class DagNode;
class MatrixReal;
template <class variableType> class StochasticNode;
    
    /**
     * The sliding operator.
     *
     * This move randomly picks an element of a matrix of positive real numbers.
     * That means, that we randomly pick the i-th row and j-th column with equal probability.
     * Then, we propose a sliding distance and slide the value.
     * The actual sliding distance is computed by delta = delta * ( u - 0.5 )
     * where u ~ Uniform(0,1).
     * The proposal is thus m[i][j] += delta * ( u - 0.5 )
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Nicolas Lartillot)
     * @since 2009-09-08, version 1.0
     *
     */
    class MatrixRealSymmetricSingleElementSlidingProposal : public Proposal {
        
    public:
        MatrixRealSymmetricSingleElementSlidingProposal( StochasticNode<MatrixReal> *n, double l);                                                                      //!<  constructor
        
        // Basic utility functions
        void                                                cleanProposal(void);                                                                //!< Clean up proposal
        MatrixRealSymmetricSingleElementSlidingProposal*    clone(void) const;                                                                  //!< Clone object
        double                                              doProposal(void);                                                                   //!< Perform proposal
        const std::string&                                  getProposalName(void) const;                                                        //!< Get the name of the proposal for summary printing
        double                                              getProposalTuningParameter(void) const;
        void                                                printParameterSummary(std::ostream &o, bool name_only) const;                                       //!< Print the parameter summary
        void                                                prepareProposal(void);                                                              //!< Prepare the proposal
        void                                                setProposalTuningParameter(double tp);
        void                                                tune(double r);                                                                     //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                                undoProposal(void);                                                                 //!< Reject the proposal
        
    protected:
        void                                                swapNodeInternal(DagNode *oldN, DagNode *newN);                                     //!< Swap the DAG nodes on which the Proposal is working on
        
    private:
        // parameters
        
        StochasticNode<MatrixReal >*                        variable;
        
        double                                              delta;                                                                             //!< The Sliding parameter of the move (larger delta -> larger proposals).
        //!< The two indices of the last modified element.
        size_t                                              indexa;
        size_t                                              indexb;
        double                                              storedValue;                                                                          //!< The value we propose.

    };
    
}

#endif

