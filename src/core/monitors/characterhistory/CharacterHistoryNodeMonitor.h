//
//  CharacterHistoryNodeMonitor.h
//  rb_mlandis
//
//  Created by Michael Landis on 10/16/13.
//  Copyright (c) 2013 Michael Landis. All rights reserved.
//

#ifndef __rb_mlandis__CharacterHistoryNodeMonitor__
#define __rb_mlandis__CharacterHistoryNodeMonitor__

#include "Monitor.h"
#include "BranchHistory.h"
#include "StochasticNode.h"
#include "TypedDagNode.h"
#include "Tree.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace RevBayesCore {
    
    class CharacterHistoryNodeMonitor : public Monitor {
        
    public:
        // Constructors and Destructors
        CharacterHistoryNodeMonitor(TypedDagNode<Tree> *t, std::vector< StochasticNode< BranchHistory >* > bh, unsigned long g, const std::string &fname, const std::string &del, bool pp=true, bool l=true, bool pr=true, bool ap=false, bool sm=true, bool sr=true);
        
        // new CharacterHistoryNodeMonitor( tau, bh_vector_stochastic, 10, filepath + "rb.tree_chars.txt", "\t"));
        
        CharacterHistoryNodeMonitor(const CharacterHistoryNodeMonitor& f);
        
        // basic methods
        CharacterHistoryNodeMonitor*        clone(void) const;                                                  //!< Clone the object
        
        // Monitor functions
        void                                monitor(unsigned long gen);                                                  //!< Monitor at generation gen
        void                                swapNode(DagNode *oldN, DagNode *newN);
        
        // FileMonitor functions
        void                                closeStream(void);                                                  //!< Close stream after finish writing
        void                                openStream(bool reopen);                                            //!< Open the stream for writing
        void                                printHeader(void);                                                  //!< Print header
        
    private:
        std::string                         buildExtendedNewick();
        std::string                         buildExtendedNewick(TopologyNode* n);
        std::string                         buildCharacterHistoryString(TopologyNode* n, std::string brEnd="child");
        
        // the stream to print
        std::fstream                        outStream;
        
        // parameters
        TypedDagNode<Tree>*                 tree;
        std::vector<StochasticNode<BranchHistory>* > branchHistories;
        std::set<DagNode *>                 nodeVariables;
        std::string                         filename;
        std::string                         separator;
        bool                                posterior;
        bool                                prior;
        bool                                likelihood;
        bool                                append;
        bool                                showMetadata;
        bool                                showRates;
        
    };
    
}


#endif /* defined(__rb_mlandis__CharacterHistoryNodeMonitor__) */
