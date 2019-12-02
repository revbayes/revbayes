/**
 * @file
 * This file contains the declaration of our tree-change event handler. Every tree has its
 * independent tree change event handler and all listener should register at this event handler
 * when they wish to be informed about tree changes.
 *
 *
 * @brief Declaration of the TreeChangeEventHandler
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

#ifndef TreeChangeEventHandler_H
#define TreeChangeEventHandler_H

#include "TopologyNode.h"

#include <set>
#include <string>

namespace RevBayesCore {
    
    class TreeChangeEventListener;
    
    class TreeChangeEventHandler {
        
    public:
        // constructors and destructors
        TreeChangeEventHandler(void);
        TreeChangeEventHandler(const TreeChangeEventHandler &h);
        virtual                                    ~TreeChangeEventHandler(void);
        
        TreeChangeEventHandler&                     operator=(const TreeChangeEventHandler &h);
        
        // public methods
        void                                        addListener(TreeChangeEventListener* l);                        //!< Add a new listener
        void                                        fire(const TopologyNode& n, const unsigned& m=0);
        const std::set<TreeChangeEventListener*>&   getListeners(void) const;
        bool                                        isListening(TreeChangeEventListener* l) const;                  //!< Is this listener listening to this tree?
        void                                        removeListener(TreeChangeEventListener* l);                     //!< Remove an existant listener
        
    private:
        std::set<TreeChangeEventListener*>          listeners;
    };

}

#endif
