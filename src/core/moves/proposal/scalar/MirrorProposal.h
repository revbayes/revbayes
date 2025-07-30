#ifndef MirrorProposal_H
#define MirrorProposal_H

#include <set>
#include <string>

#include "Proposal.h"
#include "ContinuousStochasticNode.h"

namespace RevBayesCore {

    /**
     * The sliding bactrian operator.
     *
     * A mirror draws a random number from a normal distribution centered on the opposite side of the posterior
     * using the posterior mean to find this opposite.
     * Lambda is the tuning parameter of the proposal to influence the size of the proposals, the SD of the normals.
     * See Thawornwattana et al. (2017) for more details.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Andrew Magee)
     * @since 2009-09-08, version 1.0
     *
     */
    class MirrorProposal : public Proposal {

    public:
        MirrorProposal( StochasticNode<double> *n, double l, double m0, double n0, double c0, double m, bool a, double p=0.44);                                                                      //!<  constructor

        // Basic utility functions
        void                                    cleanProposal(void);                                                                //!< Clean up proposal
        MirrorProposal*                         clone(void) const;                                                                  //!< Clone object
        double                                  doProposal(void);                                                                   //!< Perform proposal
        const std::string&                      getProposalName(void) const;                                                        //!< Get the name of the proposal for summary printing
        double                                  getProposalTuningParameter(void) const;
        void                                    printParameterSummary(std::ostream &o, bool name_only) const;                                       //!< Print the parameter summary
        void                                    prepareProposal(void);                                                              //!< Prepare the proposal
        void                                    setProposalTuningParameter(double tp);
        void                                    tune(double r);                                                                     //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                    undoProposal(void);                                                                 //!< Reject the proposal

    protected:
        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);                                     //!< Swap the DAG nodes on which the Proposal is working on

    private:
        // parameters

        StochasticNode<double>*                 variable;                                                                           //!< The variable the Proposal is working on
        double                                  storedValue;                                                                        //!< The stored value of the Proposal used for rejections.
        double                                  lambda;                                                                             //!< The value we propose.
        double                                  muHat;                                                                              //!< The currently used estimated mean of the variable from MCMC samples so far
        double                                  learnedMuHat;                                                                       //!< The estimated mean of the variable from MCMC samples so far
        double                                  sigmaHat;                                                                           //!< The currently used estimated sd of the variable from MCMC samples so far
        double                                  learnedSigmaHat;                                                                    //!< The estimated sd of the variable from MCMC samples so far

        size_t                                  waitBeforeLearning;                                                                 //!< How std::int64_t to wait before tracking empirical covariances
        size_t                                  waitBeforeUsing;                                                                    //!< How std::int64_t to wait before using the empirical covariances
        size_t                                  nTried;                                                                             //!< How many times has this move been used?
        size_t                                  updates;                                                                            //!< How many updates have been tried?
        size_t                                  maxUpdates;                                                                         //!< How many updates until we stop monitoring the covariances?

        bool                                    adaptOnly;                                                                          //!< If true, proposal variance is not tuned and user-specified value is constant. Adaptation still occurs
    };

}

#endif
