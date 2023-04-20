#ifndef EmpiricalTreeTopologyProposal_H
#define EmpiricalTreeTopologyProposal_H

#include <string>

#include "Proposal.h"
#include "StochasticNode.h"
#include "Tree.h"

namespace RevBayesCore {
    
    class TraceTree;

    /**
     * The narrow-exchange operator.
     *
     * A narrow-exchange proposal is a NNI (nearest neighbour interchange) proposal on rooted trees without changing the node age.
     * That is, we pick a random node which is not the root and neither its parent is the root.
     * Then, we try to exchange the picked node with it's uncle. This move will automatically fail if the uncle is older than the parent.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2012-07-12, version 1.0
     *
     */
    class EmpiricalTreeTopologyProposal : public Proposal {
        
    public:

        EmpiricalTreeTopologyProposal( StochasticNode<Tree> *t, const TraceTree &tt, double a );                                               //!<  constructor
        
        // Basic utility functions
        void                                    cleanProposal(void);                                        //!< Clean up proposal
        EmpiricalTreeTopologyProposal*          clone(void) const;                                          //!< Clone object
        double                                  doProposal(void);                                           //!< Perform proposal
        const std::string&                      getProposalName(void) const;                                //!< Get the name of the proposal for summary printing
        void                                    prepareProposal(void);                                      //!< Prepare the proposal
        void                                    printParameterSummary(std::ostream &o) const;               //!< Print the parameter summary
        void                                    tune(double r);                                             //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                    undoProposal(void);                                         //!< Reject the proposal
        
    protected:
        
        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);             //!< Swap the DAG nodes on which the Proposal is working on
        
        
    private:
        
        
        // member variables
        StochasticNode<Tree>*                   tree;
        std::vector<Tree>                       trees;
        
        // parameters
        double                                  alpha;
        
        // stored objects to undo proposal
        std::vector<double>                     stored_ages;
        double                                  stored_root_age;

    };
    
}

#endif

