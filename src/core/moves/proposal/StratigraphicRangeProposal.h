#ifndef StratigraphicRangeProposal_H
#define StratigraphicRangeProposal_H

#include "Proposal.h"
#include "StochasticNode.h"
#include "Tree.h"

namespace RevBayesCore {

    class FossilizedBirthDeathSpeciationProcess;

    /**
     * Resample one taxon's stratigraphic range: the first and last appearances, each drawn within
     * the bins its occurrences reported.
     *
     * Specific to dnFBDSP. A tree carries the divergence times and nothing else, so the two
     * appearances are the distribution's own state and no move on the value can reach them. The
     * matrix process keeps them in its value instead, where the generic element moves already do.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2012-07-12, version 1.0
     *
     */
    class StratigraphicRangeProposal : public Proposal {

    public:
        StratigraphicRangeProposal( StochasticNode<Tree> *n );                                      //!< Constructor

        // Basic utility functions
        bool                                    allowClamped() const override { return true; }      //!< Samples the appearances rather than the clamped tree, so it is valid on a clamped node. See #600.
        void                                    cleanProposal(void);                                //!< Clean up proposal
        StratigraphicRangeProposal*             clone(void) const;                                  //!< Clone object
        double                                  doProposal(void);                                   //!< Perform proposal
        const std::string&                      getProposalName(void) const;                        //!< Get the name of the proposal for summary printing
        double                                  getProposalTuningParameter(void) const;
        void                                    prepareProposal(void);                              //!< Prepare the proposal
        void                                    printParameterSummary(std::ostream &o, bool name_only) const; //!< Print the parameter summary
        void                                    undoProposal(void);                                 //!< Reject the proposal

    protected:

        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);     //!< Swap the DAG nodes the Proposal is working on

    private:

        FossilizedBirthDeathSpeciationProcess&  rangeProcess(void) const;                           //!< The distribution, which owns the appearances its value does not carry.

        StochasticNode<Tree>*                   variable;                                           //!< The variable the Proposal is working on
    };

}

#endif
