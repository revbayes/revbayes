#ifndef NearestNeighborInterchange_nonClockProposal_H
#define NearestNeighborInterchange_nonClockProposal_H

#include <iosfwd>

#include "Proposal.h"

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
    class NearestNeighborInterchange_nonClockProposal : public Proposal {
        
    public:
        NearestNeighborInterchange_nonClockProposal( StochasticNode<Tree> *n );                                               //!<  constructor
        
        // Basic utility functions
        void                                                cleanProposal(void);                                        //!< Clean up proposal
        NearestNeighborInterchange_nonClockProposal*        clone(void) const;                                          //!< Clone object
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
        
        
    private:
        
        // member variables
        StochasticNode<Tree>*                               tree = nullptr;
        
        // stored objects to undo proposal
        TopologyNode*                                       stored_node_A = nullptr;
        TopologyNode*                                       stored_node_B = nullptr;
        bool                                                picked_root_branch = false;
        bool                                                picked_uncle = false;
        bool                                                failed = false;
        
    };
    
}

#endif

