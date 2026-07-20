#ifndef GibbsBuddingTopologyProposal_H
#define GibbsBuddingTopologyProposal_H

#include <iosfwd>

#include "Proposal.h"
#include "StochasticNode.h"
#include "Tree.h"

namespace RevBayesCore {

    /**
     * Redraw the whole budding (asymmetric speciation) topology of a fossilized birth death speciation tree, holding the
     * ranges fixed.
     *
     * Every lineage buds off one drawn uniformly from those alive at its birth. Conditional on the
     * ranges each compatible tree carries the same density, which is what the range process states
     * as a factor of gamma per taxon, so the draw is from the exact conditional and the move is a
     * Gibbs step: the density ratio is one and the proposal is always accepted.
     *
     * The local topology moves reach these trees by rearranging one branch at a time; this one
     * lands anywhere in the compatible set in a single step.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (June Walker)
     * @since 2026-07-19, version 1.0
     *
     */
    class GibbsBuddingTopologyProposal : public Proposal {

    public:
        GibbsBuddingTopologyProposal( StochasticNode<Tree> *n );                                      //!< constructor

        // Basic utility functions
        void                                    cleanProposal(void);                                            //!< Clean up proposal
        GibbsBuddingTopologyProposal*           clone(void) const;                                              //!< Clone object
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

        StochasticNode<Tree>*                   variable;                                                       //!< The variable the Proposal is working on

        Tree                                    stored_tree;
        bool                                    failed;
    };

}

#endif
