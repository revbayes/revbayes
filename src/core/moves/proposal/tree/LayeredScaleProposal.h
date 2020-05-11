#ifndef LayeredScaleProposal_H
#define LayeredScaleProposal_H

#include <iosfwd>
#include <vector>

#include "Proposal.h"

namespace RevBayesCore {
class DagNode;
class Tree;
template <class variableType> class StochasticNode;

    /**
     * The subtree-scale operator.
     *
     * A layered scale proposal is a scaling proposal that changes the scale of all subtrees below some given age,
     * without changing the topology.
     * That is, we pick a random age between the root age and tip ages.
     * Then, we randomky pick a multiplier.
     * All subtree ages below the age are rescaled with this multiplier.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Bastien Bousau)
     * @since 2012-07-12, version 1.0
     *
     */
    class LayeredScaleProposal : public Proposal {

    public:
        LayeredScaleProposal( StochasticNode<Tree> *n, double l, double p=0.44);                            //!<  constructor

        // Basic utility functions
        void                                    cleanProposal(void);                                        //!< Clean up proposal
        LayeredScaleProposal*                   clone(void) const;                                          //!< Clone object
        double                                  doProposal(void);                                           //!< Perform proposal
        const std::string&                      getProposalName(void) const;                                //!< Get the name of the proposal for summary printing
        double                                  getProposalTuningParameter(void) const;
        void                                    prepareProposal(void);                                      //!< Prepare the proposal
        void                                    printParameterSummary(std::ostream &o, bool name_only) const;               //!< Print the parameter summary
        void                                    setProposalTuningParameter(double tp);
        void                                    tune(double r);                                             //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                    undoProposal(void);                                         //!< Reject the proposal

    protected:

        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);             //!< Swap the DAG nodes on which the Proposal is working on


    private:

        // parameters
        StochasticNode<Tree>*                   variable;                                                   //!< The variable the Proposal is working on
        double                                  lambda;                                                     //!< The scaling parameter of the Proposal

        // stored objects to undo proposal
        std::vector<double>                     storedAges;

    };

}

#endif
