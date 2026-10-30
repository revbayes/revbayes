#ifndef NarrowExchangeProposal_H
#define NarrowExchangeProposal_H

#include <iosfwd>

#include "RbVector.h"
#include "Proposal.h"

namespace RevBayesCore {

class DagNode;
class TopologyNode;
class Tree;

template <class variableType> class StochasticNode;
    
    /**
     * The narrow-exchange operator.
     *
     * A narrow-exchange proposal is a NNI (nearest neighbour interchange) on rooted trees.
     * Pick a non-root node whose parent is not the root, and exchange it with its uncle.
     * All node ages stay fixed. An ordinary uncle must be younger than the parent.
     * An exchanged sampled-ancestor tip must already match its destination parent's age.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2012-07-12, version 1.0
     *
     */
    class NarrowExchangeProposal : public Proposal {
        
    public:
        NarrowExchangeProposal( StochasticNode<Tree> *n, StochasticNode< RbVector<Tree> > *vec_s );             //!< Constructor
        
        // Basic utility functions
        void                                    cleanProposal(void);                                            //!< Clean up proposal
        NarrowExchangeProposal*                 clone(void) const;                                              //!< Clone object
        double                                  doProposal(void);                                               //!< Perform proposal
        const std::string&                      getProposalName(void) const;                                    //!< Get the name of the proposal for summary printing
        double                                  getProposalTuningParameter(void) const;
        void                                    prepareProposal(void);                                          //!< Prepare the proposal
        void                                    printParameterSummary(std::ostream &o, bool name_only) const;   //!< Print the parameter summary
        void                                    undoProposal(void);                                             //!< Reject the proposal
        
    protected:
        
        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);                 //!< Swap the DAG nodes the Proposal is working on
        
        
    private:
        
        // parameters
        StochasticNode<Tree>*                   variable;                                                       //!< The variable the proposal is working on
        StochasticNode< RbVector<Tree> >*       vector_variable;                                                //!< The alternative variable the proposal is working on

        // stored objects to undo proposal
        bool                                    failed;
        size_t                                  tree_index;
        TopologyNode*                           storedChosenNode;
        TopologyNode*                           storedUncle;
    };
    
}

#endif

