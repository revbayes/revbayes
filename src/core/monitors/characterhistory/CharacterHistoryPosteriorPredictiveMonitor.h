//
//  CharacterHistoryPosteriorPredictiveMonitor.h
//  rb_mlandis
//
//  Created by Michael Landis on 12/23/13.
//  Copyright (c) 2013 Michael Landis. All rights reserved.
//

#ifndef __rb_mlandis__CharacterHistoryPosteriorPredictiveMonitor__
#define __rb_mlandis__CharacterHistoryPosteriorPredictiveMonitor__


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
    
    class CharacterHistoryPosteriorPredictiveMonitor : public Monitor {
        
    public:
        // Constructors and Destructors
        CharacterHistoryPosteriorPredictiveMonitor(TypedDagNode<Tree> *t, std::vector< StochasticNode< BranchHistory >* > bh, std::uint64_t g, const path &fname, const std::string &del, bool pp=true, bool l=true, bool pr=true, bool ap=false, bool sm=true, bool sr=true);
        
        // new CharacterHistoryPosteriorPredictiveMonitor( tau, bh_vector_stochastic, 10, filepath + "rb.tree_chars.txt", "\t"));
        
        CharacterHistoryPosteriorPredictiveMonitor(const CharacterHistoryPosteriorPredictiveMonitor& f);
        
        // basic methods
        CharacterHistoryPosteriorPredictiveMonitor*          clone(void) const;                                                  //!< Clone the object
        
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


#endif /* defined(__rb_mlandis__CharacterHistoryPosteriorPredictiveMonitor__) */
