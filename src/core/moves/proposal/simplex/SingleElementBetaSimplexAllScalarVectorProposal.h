#ifndef SingleElementBetaSimplexAllScalarVectorProposal_H
#define SingleElementBetaSimplexAllScalarVectorProposal_H

#include <iosfwd>
#include <vector>

#include "Proposal.h"
#include "Simplex.h"

namespace RevBayesCore {
    class DagNode;
    class TopologyNode;
    class Tree;
    template <class variableType> class StochasticNode;
    
    /**
     * The subtree-scale operator.
     *
     * A subtree-scale proposal is a scaling proposal on rooted subtrees without changing the topology.
     * That is, we pick a random node which is not the root.
     * Then, we uniformly pick an age between the parent and the oldest sampled descendant.
     * The picked subtree is then scaled to this new age.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2012-07-12, version 1.0
     *
     */
    class SingleElementBetaSimplexAllScalarVectorProposal : public Proposal {
        
    public:
        SingleElementBetaSimplexAllScalarVectorProposal( StochasticNode<Simplex> *n, std::vector<StochasticNode<double> *> sv, double a, double p=0.44);                                                     //!<  constructor
        
        // Basic utility functions
        void                                                                 addIndex(size_t index);
        
        void                                                                 cleanProposal(void);                                        //!< Clean up proposal
        SingleElementBetaSimplexAllScalarVectorProposal*   clone(void) const;                                          //!< Clone object
        double                                                               doProposal(void);                                           //!< Perform proposal
        const std::string&                                                   getProposalName(void) const;                                //!< Get the name of the proposal for summary printing
        double                                                               getProposalTuningParameter(void) const;
        void                                                                 prepareProposal(void);                                      //!< Prepare the proposal
        void                                                                 printParameterSummary(std::ostream &o, bool name_only) const;               //!< Print the parameter summary
        void                                                                 setProposalTuningParameter(double tp);
        void                                                                 tune(double r);                                             //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                                                 undoProposal(void);                                         //!< Reject the proposal
        
    protected:
        
        void                                                                 swapNodeInternal(DagNode *oldN, DagNode *newN);             //!< Swap the DAG nodes on which the Proposal is working on
        
        
    private:
        
        // parameters
        StochasticNode<Simplex>*                simplex;
        std::vector<StochasticNode<double> *>   scalar_vector;
        std::set<size_t>                        indices;
        double                                  alpha;
        
        // stored objects to undo proposal
        bool                                    failed;
        bool                                    scalar_added;
        Simplex                                 stored_simplex;
        std::vector<double>                     stored_scalar_vector;
        
    };
    
}

#endif

