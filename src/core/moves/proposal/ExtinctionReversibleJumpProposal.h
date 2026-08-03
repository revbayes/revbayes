#ifndef ExtinctionReversibleJumpProposal_H
#define ExtinctionReversibleJumpProposal_H

#include "MatrixReal.h"
#include "Proposal.h"
#include "StochasticNode.h"
#include "Tree.h"

namespace RevBayesCore {

    class AbstractFossilizedBirthDeathRangeProcess;

    /**
     * Move one taxon's extinction time between the present and a time above it.
     *
     * Under rho < 1 a taxon reported extinct has a point mass at the present, where it survived
     * and went unsampled and pays 1 - rho, mixed with a density over extinction times above it.
     * The element moves are continuous, so the point has measure zero and they never propose it.
     * This jumps between the two, which is what makes rho < 1 samplable at all.
     *
     * It takes the range matrix of a dnFBDRP, or the tree of an extended dnFBDSP, where the tip
     * age is the extinction time. A non-extended tree marginalizes those out and has no point mass.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2012-07-12, version 1.0
     *
     */
    class ExtinctionReversibleJumpProposal : public Proposal {

    public:
        ExtinctionReversibleJumpProposal( StochasticNode<MatrixReal> *n );                          //!< Constructor
        ExtinctionReversibleJumpProposal( StochasticNode<Tree> *n );                                //!< Constructor, extended trees only

        // Basic utility functions
        void                                    cleanProposal(void);                                //!< Clean up proposal
        ExtinctionReversibleJumpProposal*       clone(void) const;                                  //!< Clone object
        double                                  doProposal(void);                                   //!< Perform proposal
        const std::string&                      getProposalName(void) const;                        //!< Get the name of the proposal for summary printing
        double                                  getProposalTuningParameter(void) const;
        void                                    prepareProposal(void);                              //!< Prepare the proposal
        void                                    printParameterSummary(std::ostream &o, bool name_only) const; //!< Print the parameter summary
        void                                    setProposalTuningParameter(double tp);
        void                                    tune(double r);                                     //!< Tune the proposal
        void                                    undoProposal(void);                                 //!< Reject the proposal

    protected:

        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);     //!< Swap the DAG nodes the Proposal is working on

    private:

        AbstractFossilizedBirthDeathRangeProcess& rangeProcess(void) const;                         //!< The distribution, which owns the present and the taxa.
        void                                    collectEligible(void);                              //!< Taxa the record reports extinct, fixed for the run.

        StochasticNode<MatrixReal>*             matrix = NULL;                                      //!< One of these is the variable, the other stays null.
        StochasticNode<Tree>*                   tree   = NULL;

        std::vector<size_t>                     eligible;                                           //!< Taxa reported extinct, the only ones with a point mass.
        size_t                                  stored_index = 0;                                   //!< The row touched, for the restore.
        double                                  stored_death = 0.0;
    };

}

#endif
