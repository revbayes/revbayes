#ifndef TreeScaleProposal_H
#define TreeScaleProposal_H

#include <iosfwd>
#include <vector>

#include "Proposal.h"
#include "RbVector.h"

namespace RevBayesCore {

    class DagNode;
    class Tree;
    template <class variableType> class StochasticNode;
    
    /**
     * The tree-scale operator.
     *
     * Scale internal ages together, stopping at tips and sampled-ancestor parents.
     * The Hastings factor counts only the scaled ages. An optional root-age scalar tracks
     * the tree root and is not an additional independent coordinate.
     * Alternatively, scale an independent origin age with the tree, adding one coordinate.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2012-07-12, version 1.0
     *
     */
    class TreeScaleProposal : public Proposal {
        
    public:
        TreeScaleProposal(StochasticNode<Tree>* t, StochasticNode<RbVector<Tree>>* vec_n,
                          StochasticNode<double>* r, double d, StochasticNode<double>* o = nullptr);
        
        // Basic utility functions
        void                                    cleanProposal(void);                                                //!< Clean up proposal
        TreeScaleProposal*                      clone(void) const;                                                  //!< Clone object
        double                                  doProposal(void);                                                   //!< Perform proposal
        const std::string&                      getProposalName(void) const;                                        //!< Get the name of the proposal for summary printing
        double                                  getProposalTuningParameter(void) const;
        void                                    prepareProposal(void);                                              //!< Prepare the proposal
        void                                    printParameterSummary(std::ostream &o, bool name_only) const;       //!< Print the parameter summary
        void                                    setProposalTuningParameter(double tp);
        void                                    tune(double r);                                                     //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                    undoProposal(void);                                                 //!< Reject the proposal
        
    protected:
        
        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);                     //!< Swap the DAG nodes the Proposal is working on
        
        
    private:
        
        
        // member variables
        StochasticNode<Tree>*                   tree;
        StochasticNode< RbVector<Tree> >*       vector_variable;                                                    //!< The alternative variable the proposal is working on
        StochasticNode<double>*                 rootAge;
        StochasticNode<double>*                 originAge;
        
        // parameters
        double                                  delta;
        
        // stored objects to undo proposal
        size_t                                  tree_index;
        std::vector<double>                     storedAges;
        double                                  storedRootAge;
        double                                  storedOriginAge;
    };
    
}

#endif

