#ifndef RateAgeSubtreeProposal_H
#define RateAgeSubtreeProposal_H

#include "RbVector.h"
#include "Proposal.h"
#include "StochasticNode.h"
#include "Tree.h"

#include <ostream>
#include <vector>
#include <string>

namespace RevBayesCore {
    
    /**
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @copyright GPL version 3
     * @since 2015-05-25, version 1.0
     */
    class RateAgeSubtreeProposal : public Proposal {
        
    public:
        RateAgeSubtreeProposal(StochasticNode<Tree> *n, double a, double p=0.44);                                                       //!< Constructor
        
        void                                        addRates(std::vector<StochasticNode<double> *> v);                                  //!< Add an up-scaling variable

        void                                        cleanProposal(void);                                                                //!< Clean up proposal
        RateAgeSubtreeProposal*                     clone(void) const;                                                                  //!< Clone object
        double                                      doProposal(void);                                                                   //!< Perform proposal
        const std::string&                          getProposalName(void) const;                                                        //!< Get the name of the proposal for summary printing
        double                                      getProposalTuningParameter(void) const;
        void                                        printParameterSummary(std::ostream &o, bool name_only) const;                       //!< Print the parameter summary
        void                                        prepareProposal(void);                                                              //!< Prepare the proposal
        void                                        setProposalTuningParameter(double tp);
        void                                        tune(double r);                                                                     //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                        undoProposal(void);                                                                 //!< Reject the proposal
        
    protected:
        
        void                                        swapNodeInternal(DagNode *oldN, DagNode *newN);                                     //!< Swap the DAG nodes the Proposal is working on
        
        
    private:
        
        // parameters
        StochasticNode<Tree>*                                tree;
        std::vector< std::vector<StochasticNode<double> *> > branch_rates;
        
        // stored objects to undo proposal
        TopologyNode*                                        stored_node;
        std::vector<double>                                  stored_ages;
        std::vector<double>                                  stored_branch_lengths;
        
        double                                               alpha;                                                                     //!< The scale parameter of the proposal (smaller alpha -> smaller a and b -> larger beta distr. variance -> larger proposals)
        
        
    };
    
}

#endif

