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
        void                                    cleanProposal(void);                                        //!< Clean up proposal
        MultiValueEventSlideProposal*           clone(void) const;                                          //!< Clone object
        double                                  doProposal(void);                                           //!< Perform proposal
        const std::string&                      getProposalName(void) const;                                //!< Get the name of the proposal for summary printing
        double                                  getProposalTuningParameter(void) const;
        void                                    prepareProposal(void);                                      //!< Prepare the proposal
        void                                    printParameterSummary(std::ostream &o, bool name_only) const;               //!< Print the parameter summary
        void                                    setProposalTuningParameter(double tp);
        void                                    tune(double r);                                             //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                    undoProposal(void);                                         //!< Reject the proposal
        
    protected:
        
        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);             //!< Swap the DAG nodes on which the Proposal is working on
        
    private:
        
        // parameters
        StochasticNode<MultiValueEvent>*        event_var;                                                   //!< The variable the Proposal is working on
        
        // stored objects to undo proposal
        bool                                    failed;
        std::string                             value_name;
        //        std::vector<double>                     stored_values;
        double                                  lambda;                                                                             //!< The Slide parameter of the move (larger lambda -> larger proposals).
        double                                  stored_value;
        size_t                                  stored_index;
        
    };
    
}

#endif

