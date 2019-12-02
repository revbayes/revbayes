#ifndef SpeciesTreeNodeSlideProposal_H
#define SpeciesTreeNodeSlideProposal_H

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
    class SpeciesTreeNodeSlideProposal : public Proposal {

    public:
        SpeciesTreeNodeSlideProposal(StochasticNode<Tree> *sp, StochasticNode<double> *r, double d );                   //!<  constructor
        virtual ~SpeciesTreeNodeSlideProposal(void);

        // Basic utility functions
        void                                            addGeneTree(StochasticNode<Tree> *gt);                          //!< Add a DAG Node holding a gene tree on which this move should operate on
        void                                            cleanProposal(void);                                                //!< Clean up proposal
        SpeciesTreeNodeSlideProposal*                   clone(void) const;                                                  //!< Clone object
        double                                          doProposal(void);                                                   //!< Perform proposal
        const std::string&                              getProposalName(void) const;                                        //!< Get the name of the proposal for summary printing
        double                                          getProposalTuningParameter(void) const;
        void                                            prepareProposal(void);                                              //!< Prepare the proposal
        void                                            printParameterSummary(std::ostream &o, bool name_only) const;                       //!< Print the parameter summary
        void                                            removeGeneTree(StochasticNode<Tree> *gt);                       //!< Remove a DAG Node holding a gene tree on which this move should operate on
        void                                            setProposalTuningParameter(double tp);
        void                                            tune(double r);                                                     //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                            undoProposal(void);                                                 //!< Reject the proposal

    protected:

        void                                            swapNodeInternal(DagNode *oldN, DagNode *newN);                     //!< Swap the DAG nodes on which the Proposal is working on


    private:

        // helper methods
        void                                            mauCanonical(Tree &tree, std::vector<TopologyNode*> &order, std::vector<bool>& wasSwapped);
        size_t                                          mauCanonicalSub(Tree &tree, TopologyNode *node, size_t loc, std::vector<TopologyNode*> &order, std::vector<bool>& wasSwaped);
        void                                            mauReconstruct(Tree &tree, std::vector<TopologyNode*> &order, std::vector<bool>&swapped);
        TopologyNode*                                   mauReconstructSub(Tree &tree, size_t from, size_t to, std::vector<TopologyNode*> &order, std::vector<bool>&wasSwaped);
        void                                            fillPreorderIndices(Tree &t, std::vector<size_t> &);
        size_t                                          fillPreorderIndices(TopologyNode &n, size_t loc, std::vector<size_t> &);

        // parameters
        StochasticNode<Tree>*                           speciesTree;                                                        //!< The variable the Proposal is working on
        std::vector< StochasticNode<Tree> *>            geneTrees;
        StochasticNode<double>*                         rootAge;

        // parameters
        double                                          delta;

        // stored objects to undo proposal
        double                                          storedAge;
        TopologyNode*                                   storedNode;
        std::vector<bool>                               swappedNodes;
        std::vector<TopologyNode*>                      orderedNodes;

//        std::vector<size_t>                             preOrderIndexBefore;
//        std::vector<size_t>                             preOrderIndexAfter;
    };

}

#endif
