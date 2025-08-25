#ifndef NodeTimeScaleProposal_H
#define NodeTimeScaleProposal_H

#include <string>

#include "Proposal.h"
#include "StochasticNode.h"
#include "Tree.h"

namespace RevBayesCore {
    
    class TopologyNode;
    
    /**
     * The node-age scale proposal operator using a Uniform distribution.
     *
     *
     * This node-age proposal works on rooted subtrees without changing the topology.
     * That is, we pick a random node which is not the root.
     * Then a scaling proposal draws a random uniform number u ~ unif (-0.5,0.5)
     * and scales the current node age by a scaling factor
     * sf = exp( lambda * u )
     * where lambda is the tuning parameter of the Proposal to influence the size of the proposals.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Will Freyman)
     * @since 2012-07-12, version 1.0
     *
     */
    class NodeTimeScaleProposal : public Proposal {
        
    public:
        NodeTimeScaleProposal( StochasticNode<Tree> *n, double l);                                               //!<  constructor
        
        // Basic utility functions
        void                                    cleanProposal(void);                                        //!< Clean up proposal
        NodeTimeScaleProposal*                  clone(void) const;                                          //!< Clone object
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
        double                                  lambda;
        
        // stored objects to undo proposal
        TopologyNode*                           storedNode = NULL;
        double                                  storedAge;
        
    };
    
}

#endif
