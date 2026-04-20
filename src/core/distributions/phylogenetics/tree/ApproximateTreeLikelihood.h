#ifndef ApproximateTreeLikelihood_H
#define ApproximateTreeLikelihood_H

#include <cstddef>
#include <vector>

#include "MatrixReal.h"
#include "Tree.h"
#include "TreeChangeEventListener.h"
#include "TypedDistribution.h"
#include "TypedDagNode.h"
#include "Clade.h"
#include "TopologyNode.h"

namespace RevBayesCore {
class DagNode;


    enum TRANSFORMATION { NONE, LOG, SQRT, ARCSIN };

    class ApproximateTreeLikelihood : public TypedDistribution<Tree>, public TreeChangeEventListener {

    public:
        ApproximateTreeLikelihood(const TypedDagNode<Tree>* tt, TypedDagNode< RbVector<double> > *br, const RbVector<double> *gr, const MatrixReal *h, TRANSFORMATION tr);
        virtual                                            ~ApproximateTreeLikelihood(void);                                        //!< Virtual destructor

        // public member functions
        ApproximateTreeLikelihood*                          clone(void) const;                                                      //!< Create an independent clone
        double                                              computeLnProbability(void);
        virtual void                                        fireTreeChangeEvent(const TopologyNode &n, const unsigned& m=0);        //!< This node was changed in the tree
        void                                                redrawValue(void);
        virtual void                                        setValue(Tree *v, bool f=false);                                        //!< Set the current value, e.g. attach an observation (clamp)

    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);        //!< Swap a parameter

    private:

        // helper functions
        bool                                                checkTopologyMatch(void) const;
        std::vector<double>                                 computeBranchLengths(void);
        void                                                initializeMapping(void);                                                //!< Map each time-tree branch by split to indices in mle_branch_lengths / gradients / Hessian
        void                                                simulateTree(void);
        void                                                transformBranchLengths(std::vector<double>& bl);
        
        virtual void                                        keepSpecialization(const DagNode* affecter);
        virtual void                                        restoreSpecialization(const DagNode *restorer);
        virtual void                                        touchSpecialization(const DagNode *toucher, bool touchAll);
        
        // members
        const TypedDagNode<Tree>*                           time_tree;
        const TypedDagNode<RbVector<double> >*              branch_rates;
        const RbVector<double>*                             gradients;
        const MatrixReal*                                   hessian;

        TRANSFORMATION                                      transform;
        bool                                                topology_match_checked;
        std::vector<double>                                 mle_branch_lengths;
        bool                                                mle_preorder;                                                           //!< Are the nodes of the MLE tree indexed by preorder traversal? True for setValue(), false for simulateTree()
        std::vector<size_t>                                 internal_to_mle_index;                                                  //!< A mapping vector: indices give the non-root nodes of the proposed (internal) time tree; values give the corresponding node indices in the unrooted MLE tree

    };

}

#endif
