#ifndef TreeAssemblyFunction_H
#define TreeAssemblyFunction_H

#include <cstddef>
#include <set>
#include <vector>

#include "Tree.h"
#include "TreeChangeEventListener.h"
#include "TypedFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;
    
    class TreeAssemblyFunction : public TypedFunction<Tree>, public TreeChangeEventListener {
        
    public:
        TreeAssemblyFunction(const TypedDagNode<Tree> *t, const TypedDagNode< RbVector<double> > *b);
        TreeAssemblyFunction(const TreeAssemblyFunction &f);
        virtual                                            ~TreeAssemblyFunction(void);                                                         //!< Virtual destructor
        
        // overloaded operators
        TreeAssemblyFunction&                               operator=(const TreeAssemblyFunction &d);
        
        // public member functions
        TreeAssemblyFunction*                               clone(void) const;                                                                  //!< Create an independent clone
        void                                                keep(const DagNode* affecter);
        void                                                restore(const DagNode *restorer);
        void                                                reInitialized(void);                                                                //!< The arguments have been re-initialized
        void                                                touch(const DagNode *toucher );
        void                                                update(void);
        void                                                fireTreeChangeEvent(const TopologyNode &n, const unsigned& m = 0);                   //!< The input tree says one of its nodes changed.
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                    //!< Implementation of swaping parameters
        
    private:
        
        // members
        void                                                copyTopology(void);                                                               //!< Copy the input's topology, as a tree of branch lengths.
        void                                                listenTo(const Tree *t);

        const TypedDagNode<Tree>*                           tau = nullptr;
        const TypedDagNode< RbVector<double> >*             brlen = nullptr;
        std::set<size_t>                                    touchedNodeIndices;
        bool                                                brlenFlagDirty = true;
        const Tree*                                         heard_from = NULL;                                                                  //!< The tree we listen to, which dnEmpiricalTree can replace.
        std::vector<size_t>                                 event_nodes;                                                                        //!< Nodes the input fired a tree change event for, to be fired on ours once it has caught up.
        bool                                                heard_event = false;
        bool                                                tau_touched = true;
        bool                                                shape_moved = true;                                                                 //!< Whether our copy's shape is behind the input's, until the topology is copied again.
    };
    
}

#endif
