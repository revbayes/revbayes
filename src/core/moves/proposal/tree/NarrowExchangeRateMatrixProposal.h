#ifndef NarrowExchangeRateMatrixProposal_H
#define NarrowExchangeRateMatrixProposal_H

#include <string>

#include "Proposal.h"
#include "RateGenerator.h"
#include "StochasticNode.h"
#include "Tree.h"

namespace RevBayesCore {
    
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
    class NarrowExchangeRateMatrixProposal : public Proposal {
        
    public:
        NarrowExchangeRateMatrixProposal( StochasticNode<Tree> *t, const std::vector<StochasticNode<RateGenerator> *> &rm);                                               //!<  constructor
        
        // Basic utility functions
        void                                    cleanProposal(void);                                        //!< Clean up proposal
        NarrowExchangeRateMatrixProposal*       clone(void) const;                                          //!< Clone object
        double                                  doProposal(void);                                           //!< Perform proposal
        const std::string&                      getProposalName(void) const;                                //!< Get the name of the proposal for summary printing
        void                                    prepareProposal(void);                                      //!< Prepare the proposal
        void                                    printParameterSummary(std::ostream &o, bool name_only) const;               //!< Print the parameter summary
        void                                    undoProposal(void);                                         //!< Reject the proposal
        
    protected:
        
        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);             //!< Swap the DAG nodes on which the Proposal is working on
        double                                  lnProposalProbabilityRateMatrix(const TopologyNode &n, size_t &oc, bool update);
        void                                    undoRateMatrixProposal(const TopologyNode &n, size_t i);

    private:
        
        // parameters
        StochasticNode<Tree>*                           tree;                                                   //!< The variable the Proposal is working on
        std::vector<StochasticNode<RateGenerator> *>    rate_matrices;
        
        // stored objects to undo proposal
        bool                                            failed;
        TopologyNode*                                   stored_choosen_node;
        TopologyNode*                                   stored_uncle;
        size_t                                          stored_index_parent;
        size_t                                          stored_index_node;
        size_t                                          stored_index_uncle;
        size_t                                          stored_index_brother;
        size_t                                          proposed_index_parent;
        size_t                                          proposed_index_node;
        size_t                                          proposed_index_uncle;
        size_t                                          proposed_index_brother;
    };
    
}

#endif

