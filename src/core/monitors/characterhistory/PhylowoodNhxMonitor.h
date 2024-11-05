#ifndef PhylowoodNhxMonitor_H
#define PhylowoodNhxMonitor_H

#include <cstddef>
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
    
    class PhylowoodNhxMonitor : public Monitor {
        
    public:
        // Constructors and Destructors
        PhylowoodNhxMonitor(TypedDagNode<Tree> *t, std::vector< StochasticNode< BranchHistory >* > bh, std::vector<std::vector<double> > gc, std::uint64_t g, std::uint64_t mg, int burn, const path &fname, const std::string &del, bool pp=true, bool l=true, bool pr=true, bool ap=false, bool sm=true, bool sr=true);
        
        // new PhylowoodNhxMonitor( tau, bh_vector_stochastic, 10, filepath + "rb.tree_chars.txt", "\t"));
        
        PhylowoodNhxMonitor(const PhylowoodNhxMonitor& f);
        
        // basic methods
        PhylowoodNhxMonitor*          clone(void) const;                                                  //!< Clone the object
        
        // Monitor functions
        void                                monitor(std::uint64_t gen);                                                  //!< Monitor at generation gen
        void                                swapNode(DagNode *oldN, DagNode *newN);
        
        // FileMonitor functions
        void                                closeStream(void);                                                  //!< Close stream after finish writing
        void                                openStream(bool reopen);                                            //!< Open the stream for writing
        void                                printHeader(void);                                                  //!< Print header
        std::vector<unsigned int>           getChildCharacterCounts(size_t idx);
        std::vector<unsigned int>           getParentCharacterCounts(size_t idx);
        std::int64_t                                getNumSamples(void);
        
    private:
        std::string                         buildExtendedNewick();
        std::string                         buildExtendedNewick(TopologyNode* n);
        std::string                         buildCharacterHistoryString(TopologyNode* n, std::string brEnd="child");
        void updateCharacterCounts(TopologyNode* n, std::string brEnd="child");
        std::string                         buildNhxString(void);
        
        // the stream to print
        std::fstream                        outStream;
        
        // parameters
        TypedDagNode<Tree>*                 tree;
        std::vector<StochasticNode<BranchHistory>* > branchHistories;
        std::set<DagNode *>                 nodeVariables;
        std::vector<std::vector<double> > geographicCoordinates;
        
        std::vector<std::vector<unsigned int> > parentCharacterCounts;
        std::vector<std::vector<unsigned int> > childCharacterCounts;
        
        size_t numHistories;
        size_t numCharacters;
        
        path                                filename;
        std::string                         separator;
        bool                                posterior;
        bool                                prior;
        bool                                likelihood;
        bool                                append;
        bool                                showMetadata;
        bool                                showRates;
        std::int64_t                                numSamples;
        std::uint64_t                       maxGen;
        std::int64_t                                burn;
        
    };
    
}


#endif
