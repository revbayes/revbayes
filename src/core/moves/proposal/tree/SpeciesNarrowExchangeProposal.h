#ifndef SpeciesNarrowExchangeProposal_H
#define SpeciesNarrowExchangeProposal_H

#include <string>
#include <set>

#include "Proposal.h"
#include "StochasticNode.h"
#include "Tree.h"

namespace RevBayesCore {
    
    /**
     * The species-narrow-exchange operator.
     *
     * A species-node-time-slide proposal is a sliding window proposal on rooted subtrees without changing the topology.
     * That is, we pick a random node which is not the root.
     * Then, we uniformly propose a new age centered around the current age plus-minus epsilon.
     * All gene-trees that are present in the population will be scaled too!
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2015-06-23, version 1.0
     *
     */
    class SpeciesNarrowExchangeProposal : public Proposal {
        
    public:
        SpeciesNarrowExchangeProposal(StochasticNode<Tree> *sp);                                                                    //!<  constructor
        
        // Basic utility functions
        void                                            addGeneTree(StochasticNode<Tree> *gt);                                      //!< Add a DAG Node holding a gene tree on which this move should operate on
        void                                            cleanProposal(void);                                                            //!< Clean up proposal
        SpeciesNarrowExchangeProposal*                  clone(void) const;                                                              //!< Clone object
        double                                          doProposal(void);                                                               //!< Perform proposal
        const std::string&                              getProposalName(void) const;                                                    //!< Get the name of the proposal for summary printing
        double                                          getProposalTuningParameter(void) const;
        void                                            prepareProposal(void);                                                          //!< Prepare the proposal
        void                                            printParameterSummary(std::ostream &o, bool name_only) const;                                   //!< Print the parameter summary
        void                                            removeGeneTree(StochasticNode<Tree> *gt);                                   //!< Remove a DAG Node holding a gene tree on which this move should operate on
        void                                            setProposalTuningParameter(double tp);
        void                                            tune(double r);                                                                 //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                            undoProposal(void);                                                             //!< Reject the proposal
        
    protected:
        
        std::vector<TopologyNode*>                      getNodesToChange( Tree &tau, TopologyNode &n, TopologyNode &b );
        std::set<TopologyNode*>                         getNodesInPopulation( Tree &tau, TopologyNode &n );
        std::set<TopologyNode*>                         getOldestSubtreesNodesInPopulation( Tree &tau, TopologyNode &n );
        std::set<TopologyNode*>                         getPossibleSiblings( TopologyNode *n, const std::set<TopologyNode*> cand );
        
        void                                            prune(TopologyNode *n, TopologyNode *c);
        void                                            regraft(TopologyNode *n, TopologyNode *c);
        
        void                                            swapNodeInternal(DagNode *oldN, DagNode *newN);                                 //!< Swap the DAG nodes on which the Proposal is working on
        
        
    private:
        
        // parameters
        StochasticNode<Tree>*                           speciesTree;                                                                    //!< The variable the Proposal is working on
        std::vector< StochasticNode<Tree> *>            geneTrees;
        
        // stored objects to undo proposal
        bool                                            failed;
        TopologyNode*                                   storedChoosenNode;
        TopologyNode*                                   storedUncle;
        std::vector< TopologyNode *>                    storedGeneTreeNodes;
        std::vector< TopologyNode *>                    storedOldBrothers;
        
    };
    
}

#endif

