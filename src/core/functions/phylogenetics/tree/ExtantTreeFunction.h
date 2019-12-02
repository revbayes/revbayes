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
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <map>
#include <set>

namespace RevBayesCore {
    
    class ExtantTreeFunction : public TypedFunction<Tree> {
        
    public:
        ExtantTreeFunction(const TypedDagNode<Tree> *t);
        virtual                                            ~ExtantTreeFunction(void);                                                         //!< Virtual destructor
        
        // public member functions
        ExtantTreeFunction*                                 clone(void) const;                                                                  //!< Create an independent clone
        void                                                keep(DagNode* affecter);
        void                                                restore(DagNode *restorer);
        void                                                reInitialized(void);                                                                //!< The arguments have been re-initialized
        void                                                touch(DagNode *toucher );
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
