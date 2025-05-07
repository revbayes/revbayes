#ifndef MultiValueEventSlideProposal_H
#define MultiValueEventSlideProposal_H

#include <cstddef>
#include <iosfwd>

#include "Proposal.h"

namespace RevBayesCore {
class DagNode;
class MultiValueEvent;
template <class variableType> class StochasticNode;
    
    /**
     * The node-age slide proposal operator using a Uniform distribution.
     *
     * This node-age proposal is a Uniform-sliding proposal on rooted subtrees without changing the topology.
     * That is, we pick a random node which is not the root.
     * Then, we pick an age between the parent and the oldest sampled descendant drawn from a Uniform distribution centered around the current age.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2012-07-12, version 1.0
     *
     */
    class MultiValueEventSlideProposal : public Proposal {
        
    public:
        MultiValueEventSlideProposal( StochasticNode<MultiValueEvent> *n, const std::string &vn, double l );                                                       //!<  constructor
        
        // Basic utility functions
        bool                                    allowClamped() const override { return true; }                         //!< Proposal doesn't change the tree, but changes parameters describing the process that generates the tree. See #600
        void                                    cleanProposal(void) override;                                          //!< Clean up proposal
        MultiValueEventSlideProposal*           clone(void) const override;                                            //!< Clone object
        double                                  doProposal(void) override;                                             //!< Perform proposal
        const std::string&                      getProposalName(void) const override;                                  //!< Get the name of the proposal for summary printing
        double                                  getProposalTuningParameter(void) const override;
        void                                    prepareProposal(void) override;                                        //!< Prepare the proposal
        void                                    printParameterSummary(std::ostream &o, bool name_only) const override; //!< Print the parameter summary
        void                                    setProposalTuningParameter(double tp) override;
        void                                    tune(double r) override;                                               //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                    undoProposal(void) override;                                           //!< Reject the proposal
        
    protected:
        
        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN) override;               //!< Swap the DAG nodes on which the Proposal is working on
        
    private:
        
        // parameters
        StochasticNode<MultiValueEvent>*        event_var;                                                             //!< The variable the proposal is working on
        
        // stored objects to undo proposal
        bool                                    failed;
        std::string                             value_name;
        //        std::vector<double>                     stored_values;
        double                                  lambda;                                                                //!< The slide parameter of the move (larger lambda -> larger proposals).
        double                                  stored_value;
        size_t                                  stored_index;
        
    };
    
}

#endif

