#ifndef StartingTreeSimulator_H
#define StartingTreeSimulator_H

#include <vector>
#include <set>
#include "RbConstants.h"


namespace RevBayesCore {
class Clade;
class Taxon;
class TopologyNode;
class Tree;
    
    /**
     * This class provides a starting tree simulator that conforms to some clade constraints.
     *
     * This class currently has only one functionality,
     * to simulate a tree with given clade and clade age constraints.
     */
    class StartingTreeSimulator {
        
    public:
        
//        enum SIM_CONDITION { TIME, SURVIVAL, ROOT };
        
        StartingTreeSimulator();
        
        Tree*                                   simulateTree( const std::vector<Taxon> &taxa, const std::vector<Clade> &constraints ) const;
        
    private:
        
        void                                    simulateClade( std::set<TopologyNode*> &nodes, double max_clade_age = RbConstants::Double::inf) const;
        
    };
    
}


#endif
