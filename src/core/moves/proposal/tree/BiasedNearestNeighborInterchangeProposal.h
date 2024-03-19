#ifndef BiasedNearestNeighborInterchangeProposal_H
#define BiasedNearestNeighborInterchangeProposal_H

#include <iosfwd>
#include <map>

#include "Proposal.h"
#include "RbBitSet.h"

namespace RevBayesCore {
class DagNode;
class TopologyNode;
class Tree;
template <class variableType> class StochasticNode;
    
    /**
     * The subtree-prune-and-regraft operator.
     *
     * A subtree-prune-and-regraft proposal is a SPR proposal on unrooted trees without changing the branch lengths.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2012-07-12, version 1.0
     *
     */
    class BiasedNearestNeighborInterchangeProposal : public Proposal {
        
    public:
        BiasedNearestNeighborInterchangeProposal( StochasticNode<Tree> *n );                                               //!<  constructor
        
        // Basic utility functions
        void                                                cleanProposal(void);                                        //!< Clean up proposal
        BiasedNearestNeighborInterchangeProposal*           clone(void) const;                                          //!< Clone object
        double                                              doProposal(void);                                           //!< Perform proposal
        const std::string&                                  getProposalName(void) const;                                //!< Get the name of the proposal for summary printing
        double                                              getProposalTuningParameter(void) const;
        void                                                prepareProposal(void);                                      //!< Prepare the proposal
        void                                                printParameterSummary(std::ostream &o, bool name_only) const;               //!< Print the parameter summary
        void                                                setProposalTuningParameter(double tp);
        void                                                tune(double r);                                             //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                                undoProposal(void);                                         //!< Reject the proposal
        
    protected:
        
        void                                                swapNodeInternal(DagNode *oldN, DagNode *newN);             //!< Swap the DAG nodes on which the Proposal is working on
        void                                                storeTree(void);
        
    private:
        
        // member variables
        StochasticNode<Tree>*                               tree;
        size_t                                              attempts;
        size_t                                              attempts_before_learning;
        size_t                                              attempts_before_using;
        size_t                                              attempts_to_learning;
        std::map<RbBitSet, size_t>                          clade_frequencies;
        
        // stored objects to undo proposal
        TopologyNode*                                       stored_node_A;
        TopologyNode*                                       stored_node_B;
        bool                                                picked_root_branch;
        bool                                                picked_uncle;
        
        
    };
    
}

#endif

