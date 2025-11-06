//
//  CharacterHistoryNodeMonitor.h
//  rb_mlandis
//
//  Created by Michael Landis on 10/16/13.
//  Copyright (c) 2013 Michael Landis. All rights reserved.
//

#ifndef __rb_mlandis__CharacterHistoryNodeMonitor__
#define __rb_mlandis__CharacterHistoryNodeMonitor__

#include <fstream>
#include <vector>
#include <set>

#include "Monitor.h"
#include "RbFileManager.h"

namespace RevBayesCore {
class BranchHistory;
class DagNode;
class TopologyNode;
class Tree;
template <class valueType> class TypedDagNode;
template <class variableType> class StochasticNode;
    
    class CharacterHistoryNodeMonitor : public Monitor {
        
    public:
        // Constructors and Destructors
        CharacterHistoryNodeMonitor(TypedDagNode<Tree> *t, std::vector< StochasticNode< BranchHistory >* > bh, std::uint64_t g, const path& fname, const std::string &del, bool pp=true, bool l=true, bool pr=true, bool ap=false, bool sm=true, bool sr=true);
        
        // new CharacterHistoryNodeMonitor( tau, bh_vector_stochastic, 10, filepath + "rb.tree_chars.txt", "\t"));
        
        CharacterHistoryNodeMonitor(const CharacterHistoryNodeMonitor& f);
        
        // basic methods
        CharacterHistoryNodeMonitor*        clone(void) const;                                                  //!< Clone the object
        
        // Monitor functions
        void                                monitor(std::uint64_t gen);                                                  //!< Monitor at generation gen
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
        path                                filename;
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
