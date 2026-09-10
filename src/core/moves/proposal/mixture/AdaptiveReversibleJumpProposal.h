#ifndef AdaptiveReversibleJumpProposal_H
#define AdaptiveReversibleJumpProposal_H

#include <set>
#include <string>

#include "Proposal.h"
#include "StochasticNode.h"

namespace RevBayesCore {
    
    /**
     * The adaptive reversible jump (RJ) proposal to switch between two elements of an RJ mixture.
     * The proposal updates or "learns" two of its member variables, sampled_mean and sampled_var
     * (representing the mean and variance of the normal proposal distribution), based on
     * previously collected samples. These values should stabilize as more samples are collected,
     * following the same logic as other adaptive moves (e.g., AVMVNProposal).
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2022-12-05, version 1.2
     *
     */
    class AdaptiveReversibleJumpProposal : public Proposal {
        
    public:
        enum PROPOSAL_DISTRIBUTION { NORMAL, GAMMA, LOGNORMAL };

        AdaptiveReversibleJumpProposal( StochasticNode<double> *n, size_t n0, size_t c0, size_t ue);                            //!< Constructor
        AdaptiveReversibleJumpProposal( const AdaptiveReversibleJumpProposal &p );
        virtual ~AdaptiveReversibleJumpProposal();
        
        AdaptiveReversibleJumpProposal&      operator=(const AdaptiveReversibleJumpProposal& p);
        
        // Basic utility functions
        void                                cleanProposal(void);                                                                //!< Clean up proposal
        AdaptiveReversibleJumpProposal*     clone(void) const;                                                                  //!< Clone object
        double                              doProposal(void);                                                                   //!< Perform proposal
        const std::string&                  getProposalName(void) const;                                                        //!< Get the name of the proposal for summary printing
        double                              getProposalTuningParameter(void) const;
        void                                prepareProposal(void);                                                              //!< Prepare the proposal
        void                                printParameterSummary(std::ostream &o, bool name_only) const;                       //!< Print the parameter summary
        void                                undoProposal(void);                                                                 //!< Reject the proposal
        
    protected:
        
        void                                swapNodeInternal(DagNode *oldN, DagNode *newN);                                     //!< Swap the DAG nodes the Proposal is working on
        
        
    private:
        // parameters
        
        StochasticNode<double>*             variable;                                                                           //!< The variable the Proposal is working on
        double                              stored_value;                                                                       //!< The stored value of the Proposal used for rejections.
        size_t                              stored_index;
        
        size_t                              wait_before_learning;                                                               //!< How long to wait before tracking empirical covariances
        size_t                              wait_before_using;                                                                  //!< How long to wait before using the empirical covariances
        size_t                              num_tried;                                                                          //!< How many times has this move been used?
        size_t                              updates_every;                                                                      //!< How frequently should we store the values for updating?
        std::vector<double>                 sampled_values;                                                                     //!< The sampled values used for updating the mean and variance
        double                              sampled_mean;                                                                       //!< The sampled mean
        double                              sampled_var;                                                                        //!< The sampled variance
        PROPOSAL_DISTRIBUTION               proposal_distribution;                                                              //!< The shape of the proposal distribution to use
    };
    
}


#endif

