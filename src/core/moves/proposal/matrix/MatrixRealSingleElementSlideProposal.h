#ifndef MatrixRealSingleElementSlideProposal_H
#define MatrixRealSingleElementSlideProposal_H

#include <cstddef>
#include <iosfwd>

#include "Proposal.h"

namespace RevBayesCore {
class DagNode;
class MatrixReal;
template <class valueType> class RbVector;
template <class variableType> class StochasticNode;
    
    /**
     * The sliding operator.
     *
     * This move randomly picks an element of a matrix of positive real numbers.
     * That means, that we randomly pick the i-th row and j-th column with equal probability.
     * Then, we propose a sliding distance and slide the value.
     * The actual sliding distance is computed by delta = lambda * ( u - 0.5 )
     * where u ~ Uniform(0,1).
     * The proposal is thus m[i][j] += lambda * ( u - 0.5 )
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Nicolas Lartillot)
     * @since 2009-09-08, version 1.0
     *
     */
    class MatrixRealSingleElementSlideProposal : public Proposal {
        
    public:
        void                                    setRow(size_t r) { restrict_row = r; }      //!< Confine the move to one row
        void                                    setColumn(size_t c) { restrict_col = c; }   //!< Confine the move to one column
        MatrixRealSingleElementSlideProposal( StochasticNode<MatrixReal> *n, double l, bool s = false);                                                                      //!< Constructor
        MatrixRealSingleElementSlideProposal( StochasticNode<RbVector<RbVector<double> > > *n, double l, bool s = false);
        
        // Basic utility functions
        void                                    cleanProposal(void);                                                                //!< Clean up proposal
        MatrixRealSingleElementSlideProposal*   clone(void) const;                                                                  //!< Clone object
        double                                  doProposal(void);                                                                   //!< Perform proposal
        const std::string&                      getProposalName(void) const;                                                        //!< Get the name of the proposal for summary printing
        double                                  getProposalTuningParameter(void) const;
        void                                    printParameterSummary(std::ostream &o, bool name_only) const;                       //!< Print the parameter summary
        void                                    prepareProposal(void);                                                              //!< Prepare the proposal
        void                                    setProposalTuningParameter(double tp);
        void                                    tune(double r);                                                                     //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                    undoProposal(void);                                                                 //!< Reject the proposal
        
    protected:
        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);                                     //!< Swap the DAG nodes the Proposal is working on
        
    private:
        // parameters
        
        StochasticNode<RbVector<RbVector<double> > >* array;
        StochasticNode<MatrixReal>*                   matrix;
        
        double                                  delta;                                                                             //!< The Sliding parameter of the move (larger delta -> larger proposals).
        //!< The two indices of the last modified element.
        size_t                                  restrict_row;                                   //!< 1-based row to confine the move to; 0 draws from the whole matrix
        size_t                                  restrict_col;                                   //!< 1-based column to confine the move to; 0 draws from the whole matrix
        size_t                                  indexa;
        size_t                                  indexb;
        double                                  storedValue;                                                                       //!< The value we propose.
        bool                                    symmetric;
    };
    
}

#endif

