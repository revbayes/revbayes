#ifndef RateAgeProposal_H
#define RateAgeProposal_H

#include "RbVector.h"
#include "Proposal.h"
#include "StochasticNode.h"
#include "Tree.h"

#include <ostream>
#include <vector>
#include <string>

namespace RevBayesCore {
    
    /**
     * @brief Up-Down Scaling Proposal for several parameters jointly.
     *
     *
     * This proposal randomly scales all a set of parameter up while the other set of parameters
     * is scaled down by the same value. This should hopefully improve mixing in many cases.
     * The actual scaling factor is computed by sf = exp( lambda * ( u - 0.5 ) )
     * where u ~ Uniform(0,1).
     *
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @copyright GPL version 3
     * @since 2015-05-25, version 1.0
     *
     */
    class RateAgeProposal : public Proposal {
        
    public:
        RateAgeProposal(StochasticNode<Tree> *n, double a, double p=0.44);                                 //!< Constructor
        
        void                                        addRates(std::vector<StochasticNode<double> *> v);                         //!< Add an up-scaling variable

        void                                        cleanProposal(void);                                                                //!< Clean up proposal
        RateAgeProposal*                            clone(void) const;                                                                  //!< Clone object
        double                                      doProposal(void);                                                                   //!< Perform proposal
        const std::string&                          getProposalName(void) const;                                                        //!< Get the name of the proposal for summary printing
        double                                      getProposalTuningParameter(void) const;
        void                                        printParameterSummary(std::ostream &o, bool name_only) const;                                       //!< Print the parameter summary
        void                                        prepareProposal(void);                                                              //!< Prepare the proposal
        void                                        setProposalTuningParameter(double tp);
        void                                        tune(double r);                                                                     //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                        undoProposal(void);                                                                 //!< Reject the proposal
        
    protected:
        
        void                                        swapNodeInternal(DagNode *oldN, DagNode *newN);                                     //!< Swap the DAG nodes on which the Proposal is working on
        
        
    private:
        
        // parameters
        StochasticNode<Tree>*                                tree;
        std::vector< std::vector<StochasticNode<double> *> > branch_rates;
        
        // stored objects to undo proposal
        TopologyNode*                                        stored_node;
        double                                               stored_age;
        std::vector< RbVector<double> >                      stored_rates;

        double                                               alpha;                                                                     //!< The scale parameter of the Proposal (larger lambda -> larger proposals).
        
        
    };
    
}

#endif

