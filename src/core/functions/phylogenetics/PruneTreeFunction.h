#ifndef PruneTreeFunction_h
#define PruneTreeFunction_h

#include <vector>

#include "Tree.h"
#include "TypedFunction.h"
#include "RbBitSet.h"
#include "Taxon.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    /**
     * @brief The Prune Tree Function
     *
     * A core function for pruning taxa from a tree.
     * The class also has the following member variables:
     *  @param prune_map An RbBitSet map for the Tree::pruneTaxa() function
     *  @param tau The tree to prune
     *  @param taxa The list of taxon names to remove
     *  @param retain Boolean: if true, then retain the tips in taxa instead of removing them
     *  @param prune_fossils Boolean: if true, then prune all of the non-extant tips/sampled ancestors from the tree
     */

    class PruneTreeFunction : public TypedFunction<Tree> {
        
    public:
        PruneTreeFunction(const TypedDagNode<Tree> *t, std::vector<Taxon> taxa, bool r, bool pf);
        virtual                                            ~PruneTreeFunction(void);                                                         //!< Virtual destructor
        
        // public member functions
        PruneTreeFunction*                                  clone(void) const;                                                                  //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                    //!< Implementation of swaping parameters
        
    private:
        
        // members
        RbBitSet                                            prune_map;

        // parameters
        const TypedDagNode<Tree>*                           tau;
        std::vector<Taxon>                                  taxa;
        bool                                                retain;
        bool                                                prune_fossils;
    };
    
}

#endif /* PruneTreeFunction_h */
