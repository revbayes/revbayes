#ifndef TreeAssemblyFunction_H
#define TreeAssemblyFunction_H

#include "RbVector.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <vector>

namespace RevBayesCore {
    
    class TreeAssemblyFunction : public TypedFunction<Tree> {
        
    public:
        TreeAssemblyFunction(const TypedDagNode<Tree> *t, const TypedDagNode< RbVector<double> > *b);
        TreeAssemblyFunction(const TreeAssemblyFunction &f);
        virtual                                            ~TreeAssemblyFunction(void);                                                         //!< Virtual destructor
        
        // overloaded operators
        TreeAssemblyFunction&                               operator=(const TreeAssemblyFunction &d);
        
        // public member functions
        TreeAssemblyFunction*                               clone(void) const;                                                                  //!< Create an independent clone
        void                                                keep(DagNode* affecter);
        void                                                restore(DagNode *restorer);   
        void                                                reInitialized(void);                                                                //!< The arguments have been re-initialized
        void                                                touch(DagNode *toucher );
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                    //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode<Tree>*                           tau;
        const TypedDagNode< RbVector<double> >*             brlen;
        std::set<size_t>                                    touchedNodeIndices;
    };
    
}

#endif
