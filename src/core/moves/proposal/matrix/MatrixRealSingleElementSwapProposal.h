#ifndef MatrixRealSingleElementSwapProposal_H
#define MatrixRealSingleElementSwapProposal_H

#include <cstddef>
#include <iosfwd>

#include "Proposal.h"
#include "StochasticNode.h"

namespace RevBayesCore {

class MatrixReal;

    /**
     * Swap two elements of a matrix, either anywhere or drawn from the same row or column.
     *
     * Confining the swap to one line exchanges the same quantity between two positions, which
     * keeps the proposal meaningful when the other margin holds different things. The swap is a
     * permutation of the value, so it is symmetric and carries no Hastings ratio.
     *
     * The margin follows R's MARGIN: 1 swaps within a row, 2 within a column, 0 anywhere.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (June Walker)
     * @since 2026-07-18, version 1.0
     *
     */
    class MatrixRealSingleElementSwapProposal : public Proposal {

    public:
        MatrixRealSingleElementSwapProposal( StochasticNode<MatrixReal> *n, std::int64_t m, std::int64_t i);   //!<  constructor

        // Basic utility functions
        void                                    cleanProposal(void);                                            //!< Clean up proposal
        MatrixRealSingleElementSwapProposal*    clone(void) const;                                              //!< Clone object
        double                                  doProposal(void);                                               //!< Perform proposal
        const std::string&                      getProposalName(void) const;                                    //!< Get the name of the proposal for summary printing
        double                                  getProposalTuningParameter(void) const;
        void                                    prepareProposal(void);                                          //!< Prepare the proposal
        void                                    printParameterSummary(std::ostream &o, bool name_only) const;   //!< Print the parameter summary
        void                                    setProposalTuningParameter(double tp);
        void                                    tune(double r);                                                 //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                    undoProposal(void);                                             //!< Reject the proposal

    protected:

        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);                 //!< Swap the DAG nodes on which the Proposal is working on

    private:

        StochasticNode<MatrixReal>*             matrix;                                                         //!< The variable the Proposal is working on
        std::int64_t                            margin;                                                         //!< 0 swaps anywhere, 1 within a row, 2 within a column (R's MARGIN)
        std::int64_t                            index;                                                          //!< Line to swap within, or -1 to draw one

        size_t                                  row_a;
        size_t                                  col_a;
        size_t                                  row_b;
        size_t                                  col_b;
        bool                                    failed;
    };

}

#endif
