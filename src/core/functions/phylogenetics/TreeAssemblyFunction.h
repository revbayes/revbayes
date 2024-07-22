#ifndef TreeAssemblyFunction_H
#define TreeAssemblyFunction_H

#include <cstddef>
#include <set>

#include "Tree.h"
#include "TypedFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;
    
    class TreeAssemblyFunction : public TypedFunction<Tree> {
        
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
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                    //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode<Tree>*                           tau = nullptr;
        const TypedDagNode< RbVector<double> >*             brlen = nullptr;
        std::set<size_t>                                    touchedNodeIndices;
        // Initialize this to false to ensure that we set branch lengths on the first update()
        bool                                                touchedTopology = false;
    };
    
}

#endif
