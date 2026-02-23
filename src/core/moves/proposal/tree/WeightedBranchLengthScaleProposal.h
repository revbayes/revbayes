#ifndef WeightedBranchLengthScaleProposal_H
#define WeightedBranchLengthScaleProposal_H

#include <stddef.h>
#include <iosfwd>

#include "Proposal.h"

namespace RevBayesCore {
class DagNode;
class Tree;
template <class variableType> class StochasticNode;
    
    /**
     * The weighted branch length scale operator.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2022-08-24, version 1.4.0-preview
     *
     */
    class WeightedBranchLengthScaleProposal : public Proposal {
        
    public:
        WeightedBranchLengthScaleProposal( StochasticNode<Tree> *t, size_t n, double a );                     //!< Constructor
        
        // Basic utility functions
        void                                    cleanProposal(void);                                          //!< Clean up proposal
        WeightedBranchLengthScaleProposal*      clone(void) const;                                            //!< Clone object
        double                                  doProposal(void);                                             //!< Perform proposal
        const std::string&                      getProposalName(void) const;                                  //!< Get the name of the proposal for summary printing
        double                                  getProposalTuningParameter(void) const;
        void                                    prepareProposal(void);                                        //!< Prepare the proposal
        void                                    printParameterSummary(std::ostream &o, bool name_only) const; //!< Print the parameter summary
        void                                    undoProposal(void);                                           //!< Reject the proposal
        
    protected:
        
        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);               //!< Swap the DAG nodes the Proposal is working on
        
        
    private:
        
        
        // member variables
        StochasticNode<Tree>*                   tree;
        
        // parameters
        size_t                                  num_breaks;
        double                                  alpha;                                                        //!< Not a tuning parameter!
        std::vector<double>                     interval;

        // stored objects to undo proposal
        double                                  stored_value;
        size_t                                  stored_branch_index;
    };
    
}

#endif


