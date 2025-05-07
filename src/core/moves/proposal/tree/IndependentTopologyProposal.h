#ifndef ConstrainedEmpiricalTreeProposal_H
#define ConstrainedEmpiricalTreeProposal_H

#include <cstddef>
#include <iosfwd>
#include <vector>

#include "Proposal.h"
#include "Tree.h"
#include "Clade.h"
#include "TopologyNode.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class RbVector;
template <class variableType> class StochasticNode;
template <class variableType> class TypedDistribution;
    
    /**
     *
     * Propose a new time tree topology from an unconstrained proposal distribution
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Will Pett)
     * @since 2012-07-12, version 1.0
     *
     */
    class IndependentTopologyProposal : public Proposal {
        
    public:
        IndependentTopologyProposal( StochasticNode<Tree> *n, TypedDistribution<Tree>* d, StochasticNode< RbVector<double> > *r = NULL, Clade og = Clade() );                                               //!<  constructor
        IndependentTopologyProposal( StochasticNode<Tree> *n, TypedDistribution<Tree>* d, std::vector<StochasticNode<double>* > r, Clade og = Clade() );                                               //!<  constructor

        
        // Basic utility functions
        void                                    cleanProposal(void);                                             //!< Clean up proposal
        IndependentTopologyProposal*            clone(void) const;                                               //!< Clone object
        double                                  doProposal(void);                                                //!< Perform proposal
        const std::string&                      getProposalName(void) const;                                     //!< Get the name of the proposal for summary printing
        double                                  getProposalTuningParameter(void) const;
        void                                    prepareProposal(void);                                           //!< Prepare the proposal
        void                                    printParameterSummary(std::ostream &o, bool name_only) const;                    //!< Print the parameter summary
        void                                    setProposalTuningParameter(double tp);
        void                                    tune(double r);                                                  //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                    undoProposal(void);                                              //!< Reject the proposal
        
    protected:
        
        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);                  //!< Swap the DAG nodes on which the Proposal is working on
        std::vector<size_t>                     recursivelyRank( const TopologyNode& n, double& ln_num_rankings );  //!< Recursively choose a random subtree node ranking
        
    private:
        
        Clade                                   outgroup;
        TypedDistribution<Tree>*                proposal_distribution;
        StochasticNode< RbVector<double> >*     rates_node;
        std::vector<StochasticNode<double>* >   rates_vector;
        StochasticNode<Tree>*                   variable;

        double                                  proposal_ln_num_rankings;
        double                                  stored_ln_num_rankings;
        std::vector<double>                     stored_rates;
        Tree                                    stored_tree;
        
        bool                                    init;
    };
    
}


#endif

