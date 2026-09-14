#ifndef FossilSiteTimeSlideUniformProposal_H
#define FossilSiteTimeSlideUniformProposal_H

#include <string>

#include "Proposal.h"
#include "StochasticNode.h"
#include "Tree.h"
#include "Taxon.h"
#include "RbVector.h"
#include "RbVectorImpl.h"

namespace RevBayesCore {
    
    /**
     * This Uniform-sliding proposal changes the age of multiple tips of the tree at once (tips originating from within one fossil site), without modifying the tree topology.
     * It ensures that after this move the included tips have the same age.
     * This move is adapted from FossilTimeSlideUniformProposal.
     * For the within-site ages to be equal in the MCC and MAP summary trees, set conditionalAges to TRUE in mccTree and mapTree.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2012-07-12, version 1.0
     *
     */
    class FossilSiteTimeSlideUniformProposal : public Proposal {
        
    public:
        FossilSiteTimeSlideUniformProposal( StochasticNode<Tree> *n,
                                           TypedDagNode<double>* o,
                                           TypedDagNode<double>* ma,
                                           TypedDagNode<double>* mi,
                                           const RbVector<Taxon> &t,
                                           double l,
                                           double r=0.44);                                                      //!<  constructor
        
        // Basic utility functions
        void                                    cleanProposal(void);                                            //!< Clean up proposal
        FossilSiteTimeSlideUniformProposal*      clone(void) const;                                              //!< Clone object
        double                                  doProposal(void);                                               //!< Perform proposal
        const std::string&                      getProposalName(void) const;                                    //!< Get the name of the proposal for summary printing
        double                                  getProposalTuningParameter(void) const;
        void                                    prepareProposal(void);                                          //!< Prepare the proposal
        void                                    printParameterSummary(std::ostream &o, bool name_only) const;   //!< Print the parameter summary
        void                                    setProposalTuningParameter(double tp);
        void                                    tune(double r);                                                 //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                    undoProposal(void);                                             //!< Reject the proposal

    protected:
        
        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);                 //!< Swap the DAG nodes on which the Proposal is working on
        
        
    private:
        
        // parameters
        StochasticNode<Tree>*                   tree;                                                           //!< The variable the Proposal is working on
        TypedDagNode<double>*                   origin;
        TypedDagNode<double>*                   max;
        TypedDagNode<double>*                   min;

        size_t                                  node_index;

        // stored objects to undo proposal
        std::vector<double> stored_ages;
        std::vector<int> stored_ages_indices;

        const RbVector<Taxon>                   taxa;
        double                                  stored_age;
        double                                  lambda;                                                                             //!< The value we propose.

        bool                                    failed;
    };
    
}

#endif
