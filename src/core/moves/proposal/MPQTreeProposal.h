#ifndef MPQTreeProposal_H
#define MPQTreeProposal_H

#include <cstddef>
#include <iosfwd>
#include <gmpxx.h>

#include "Polyhedron.h"
#include "Proposal.h"
#include "RateGenerator.h"
#include "RateMatrix_MPQ.h"

namespace RevBayesCore {
class DagNode;
class Tree;
template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;
template <class variableType> class StochasticNode;
    
    /**
     * The time-reversible and non-reversible rate matrix proposal.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (John & Sebastian)
     * @since 2009-09-08, version 1.0
     *
     */
    class MPQTreeProposal : public Proposal {
        
    public:

        enum MOVE_TYPE { BRANCH_LENGTH, TREE_LENGTH, ROOT_POSITION };

        MPQTreeProposal( TypedDagNode<RateGenerator> *q, StochasticNode<Tree>* t );                                                             //!<  constructor
                                                            MPQTreeProposal( const MPQTreeProposal& p);                                         //!<  copy constructor

        // Basic utility functions
        void                                                cleanProposal(void);                                                                //!< Clean up proposal
        MPQTreeProposal*                                    clone(void) const;                                                                  //!< Clone object
        double                                              doProposal(void);                                                                   //!< Perform proposal
        const std::string&                                  getProposalName(void) const;                                                        //!< Get the name of the proposal for summary printing
        double                                              getProposalTuningParameter(void) const;
        void                                                printParameterSummary(std::ostream &o, bool name_only) const;                                       //!< Print the parameter summary
        void                                                prepareProposal(void);                                                              //!< Prepare the proposal
        void                                                setProposalTuningParameter(double tp);
        void                                                tune(double r);                                                                     //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                                undoProposal(void);                                                                 //!< Reject the proposal
        
    protected:
        void                                                swapNodeInternal(DagNode *oldN, DagNode *newN);                                     //!< Swap the DAG nodes on which the Proposal is working on
        
    private:
        double                                              updateBranchLengths(void);          //!< Update single branch
        double                                              updateTreeLength(void);
        double                                              updateRootPosition(void);
                
        // parameters
        TypedDagNode<RateGenerator>*                        q_matrix;
        StochasticNode<Tree>*                               tree;

        // tuning parameters
        double                                              tuning_branch_length;
        double                                              tuning_tree_length;
//        double                                              rev_alpha_pi;                                                                       //!< The Sliding parameter of the move (larger lambda -> larger proposals).
//        double                                              rev_alpha_er;                                                                       //!< The Sliding parameter of the move (larger lambda -> larger proposals).
//        double                                              non_rev_alpha;                                                                      //!< The Sliding parameter of the move (larger lambda -> larger proposals).
//        double                                              rj_alpha;                                                                           //!< The Sliding parameter of the move (larger lambda -> larger proposals).

        //!< The two indices of the last modified element.
//        RateMatrix_MPQ                                      stored_Q;
        double                                              stored_branch_length;
        size_t                                              stored_branch_index;
        double                                              stored_scaling_factor;
        double                                              stored_root_index;
        MOVE_TYPE                                           last_move;

//        std::vector<mpq_class>                              W;

    };
    
}

#endif

