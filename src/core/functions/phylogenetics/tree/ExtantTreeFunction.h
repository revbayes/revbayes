//
//  ExtantTreeFunction.h
//  revbayes-proj
//
//  Created by Michael Landis on 3/2/16.
//  Copyright © 2016 Michael Landis. All rights reserved.
//

#ifndef ExtantTreeFunction_h
#define ExtantTreeFunction_h

#include "Tree.h"
#include "TypedFunction.h"

namespace RevBayesCore {
class DagNode;
class TopologyNode;
template <class valueType> class TypedDagNode;
    
    class ExtantTreeFunction : public TypedFunction<Tree> {
        
    public:
        ExtantTreeFunction(const TypedDagNode<Tree> *t);
        virtual                                            ~ExtantTreeFunction(void);                                                         //!< Virtual destructor
        
        // public member functions
        ExtantTreeFunction*                                 clone(void) const;                                                                  //!< Create an independent clone
        void                                                reInitialized(void);                                                                //!< The arguments have been re-initialized
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                    //!< Implementation of swaping parameters
        
    private:

        void                                                pruneFossilsRecursively(TopologyNode* node);
        
        // members
        
        // parameters
        const TypedDagNode<Tree>*                           tau;
        
    };
    
}

#endif /* ExtantTreeFunction_h */
