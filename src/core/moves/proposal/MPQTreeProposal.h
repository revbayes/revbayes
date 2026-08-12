#ifndef MPQTreeProposal_H
#define MPQTreeProposal_H

#include <cstddef>
#include <iosfwd>
#include <gmpxx.h>
#include <set>
#include <string>

#include "Polyhedron.h"
#include "Proposal.h"
#include "RateGenerator.h"
#include "RateMatrix_MPQ.h"

namespace RevBayesCore {
class DagNode;
class Tree;
class TopologyNode;

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

        MPQTreeProposal( StochasticNode<Tree>* t, bool ur );                                                                                 //!<  constructor
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
        void                                                setVerifyRootMove(bool tf) { verify_root_move = tf; }                               //!< Check every root move by round trip
        bool                                                lastMoveWasRoot(void) const { return last_move == ROOT_POSITION; }
        
    protected:
        void                                                swapNodeInternal(DagNode *oldN, DagNode *newN);                                     //!< Swap the DAG nodes on which the Proposal is working on
        
    private:
        double                                              updateBranchLengths(void);          //!< Update single branch
        double                                              updateTreeLength(void);
        double                                              updateRootPosition(void);
        
        void                                                markNodes( std::vector<TopologyNode*>& markedNodes, TopologyNode* curr_node );
        /* Only the tree. This move used to hold the rate matrix as well, but it
           neither reads nor modifies it, and holding a node has two consequences
           that are easy to miss.

           A node must be registered with addNode for AbstractMove::swapNode to
           ever reach it, because swapNode is dispatched through the move list the
           node itself keeps. A pointer held but not registered is never swapped,
           so once the model is cloned -- which happens for every MC3 chain and
           every nruns replicate -- it still refers to the original model.

           And registering a node is not free: MetropolisHastingsMove touches every
           registered node, so registering the rate matrix would dirty the whole
           CTMC and turn each branch-length proposal into a full-tree likelihood
           recomputation rather than one path to the root.

           Holding no pointer at all avoids both. */
        StochasticNode<Tree>*                               tree;

        /* Whether this move is allowed to move the root.

           This has to be told to us; we cannot work it out for ourselves. The
           natural test would be to propose a new root and let the tree prior
           reject it, but RevBayes issue #157 (Hoehna, March 2021) reports that
           the outgroup argument of dnUniformTopology and
           dnUniformTopologyBranchLength is enforced only when the starting tree
           is built and never again during the MCMC. If that is still true, a
           root move would silently walk the root away from the outgroup and
           nothing would object. So the Rev layer has to pass this in, and the
           caller is responsible for setting it to false whenever an outgroup
           has been assigned. */
        bool                                                update_root;

        // tuning parameters
        double                                              tuning_branch_length;
        double                                              tuning_tree_length;
//        double                                              rev_alpha_pi;                                                                       //!< The Sliding parameter of the move (larger lambda -> larger proposals).
//        double                                              rev_alpha_er;                                                                       //!< The Sliding parameter of the move (larger lambda -> larger proposals).
//        double                                              non_rev_alpha;                                                                      //!< The Sliding parameter of the move (larger lambda -> larger proposals).
//        double                                              rj_alpha;                                                                           //!< The Sliding parameter of the move (larger lambda -> larger proposals).

        //!< The two indices of the last modified element.
//        RateMatrix_MPQ                                      stored_Q;
        // store variables for the root position update
        double                                              stored_first_root_branch_length;
        double                                              stored_second_root_branch_length;
//        double                                              stored_root_branch_length_fraction;
//        double                                              stored_new_root_branch_length;
        TopologyNode*                                       stored_root_node;
        
        double                                              stored_branch_length;
        size_t                                              stored_branch_index;
        double                                              stored_scaling_factor;
        size_t                                              stored_root_index;
        MOVE_TYPE                                           last_move;

        /* Set to true to have every root-position move verified by round trip:
           the Newick string is captured before the move and compared against the
           tree undoProposal hands back. The root move is the only one here that
           rearranges topology, its undo is a hand-written reversal of that
           rearrangement, and a mistake there would corrupt the tree silently
           rather than crash. Costs a string comparison per rejected root move,
           so leave it off for production runs and on for the first few thousand
           iterations of anything new. */
        bool                                                verify_root_move;
        std::string                                         stored_newick;

//        std::vector<mpq_class>                              W;

    };
    
}

#endif

