#ifndef SpeciesSubtreeScaleProposal_H
#define SpeciesSubtreeScaleProposal_H

#include <string>

#include "Proposal.h"
#include "StochasticNode.h"
#include "Tree.h"

namespace RevBayesCore {
    
    /**
     * The species-subtree-scale operator.
     *
     * A subtree-scale proposal is a scaling proposal on rooted subtrees without changing the topology.
     * That is, we pick a random node which is not the root.
     * Then, we uniformly pick an age between the parent and the oldest sampled descendant.
     * The picked subtree is then scaled to this new age.
     * All gene-trees that are present in the population will be scaled too!
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2015-06-16, version 1.0
     *
     */
    class SpeciesSubtreeScaleProposal : public Proposal {
        
    public:
        SpeciesSubtreeScaleProposal(StochasticNode<Tree> *sp);                                               //!<  constructor
        virtual ~SpeciesSubtreeScaleProposal(void);
        
        // Basic utility functions
        void                                            addGeneTree(StochasticNode<Tree> *gt);                                      //!< Add a DAG Node holding a gene tree on which this move should operate on
        void                                            cleanProposal(void);                                        //!< Clean up proposal
        SpeciesSubtreeScaleProposal*                    clone(void) const;                                          //!< Clone object
        double                                          doProposal(void);                                           //!< Perform proposal
        const std::string&                              getProposalName(void) const;                                //!< Get the name of the proposal for summary printing
        double                                          getProposalTuningParameter(void) const;
        void                                            prepareProposal(void);                                      //!< Prepare the proposal
        void                                            printParameterSummary(std::ostream &o, bool name_only) const;               //!< Print the parameter summary
        void                                            removeGeneTree(StochasticNode<Tree> *gt);                                   //!< Remove a DAG Node holding a gene tree on which this move should operate on
        void                                            setProposalTuningParameter(double tp);
        void                                            tune(double r);                                             //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                            undoProposal(void);                                         //!< Reject the proposal
        
    protected:

        std::vector<TopologyNode*>                      getOldestNodesInPopulation( Tree &tau, TopologyNode &n );
        void                                            swapNodeInternal(DagNode *oldN, DagNode *newN);             //!< Swap the DAG nodes on which the Proposal is working on
        
        
    private:
        
        // parameters
        StochasticNode<Tree>*                           speciesTree;                                                   //!< The variable the Proposal is working on
        std::vector< StochasticNode<Tree> *>            geneTrees;
       
        // stored objects to undo proposal
        TopologyNode*                                   storedNode;
        std::vector<std::vector<double> >               storedGeneTreeAges;
        std::vector<double>                             storedSpeciesTreeAges;
        
    };
    
}

#endif

