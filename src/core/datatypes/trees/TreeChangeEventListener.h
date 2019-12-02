/**
 * @file
 * This file contains the declaration of our tree-change event listener interface. 
 * Every tree has its independent tree change event handler and all listener should 
 * register at this event handler when they wish to be informed about tree changes.
 * Possible tree change-event listeners must derive/implement from this interface.
 *
 *
 * @brief Declaration of the TreeChangeEventListener
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-07-05 16:47:08 +0200 (Thu, 05 Jul 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-08-27, version 1.0
 *
 * $Id: Tree.h 1651 2012-07-05 14:47:08Z hoehna $
 */

#ifndef TreeChangeEventListener_H
#define TreeChangeEventListener_H

#include "TopologyNode.h"
#include "TreeChangeEventMessage.h"

namespace RevBayesCore {
    
//    class TreeChangeEventHandler;
    class TreeChangeEventListener {
        
    public:
        // constructors and destructors
        TreeChangeEventListener(void) {} 
        virtual                                    ~TreeChangeEventListener(void) {}
        
        // public pure virtual methods
        virtual void                                fireTreeChangeEvent(const TopologyNode &n, const unsigned& m=0) = 0;                                 //!< This node was changed in the tree
        
    };
    
}

#endif
