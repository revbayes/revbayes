#ifndef StartingTreeSimulator_H
#define StartingTreeSimulator_H

#include <vector>
#include <set>

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
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2013-04-15, version 1.0
     */
    class StartingTreeSimulator {
        
    public:
        
//        enum SIM_CONDITION { TIME, SURVIVAL, ROOT };
        
        StartingTreeSimulator();
        
        Tree*                                   simulateTree( const std::vector<Taxon> &taxa, const std::vector<Clade> &constraints ) const;
        
    private:
        
        void                                    simulateClade( std::set<TopologyNode*> &nodes) const;
        
    };
    
}


#endif
