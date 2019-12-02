#ifndef TreeScaleFunction_H
#define TreeScaleFunction_H

#include "RbVector.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <vector>

namespace RevBayesCore {
    
    class TreeScaleFunction : public TypedFunction<Tree> {
        
    public:
        TreeScaleFunction(const TypedDagNode<Tree> *t, const TypedDagNode<double> *b, std::vector<double> m);
        virtual                                            ~TreeScaleFunction(void);                                                         //!< Virtual destructor
        
        // public member functions
        TreeScaleFunction*                                  clone(void) const;                                                                  //!< Create an independent clone
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
        const TypedDagNode<double>*                         scale;
        std::vector<double>                                 tipAges;
        double                                              scaleLimit;
    };
    
}

#endif
