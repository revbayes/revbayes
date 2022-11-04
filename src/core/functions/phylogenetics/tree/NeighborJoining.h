#ifndef NeighborJoining_hpp
#define NeighborJoining_hpp

#include <stdio.h>
#include <vector>
#include <set>

namespace RevBayesCore {
class Clade;
class DistanceMatrix;
class Taxon;
class TopologyNode;
class Tree;
    
    /**
     * This class provides the NeighborJoining algorithm to construct a tree from a distance matrix.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2013-04-15, version 1.0
     */
    class NeighborJoining {
        
    public:
        
        NeighborJoining();
        
        Tree*                   constructTree( const DistanceMatrix& d ) const;
        
    private:
        
        TopologyNode*           constructTreeRecursively( std::vector<TopologyNode*> &nodes, MatrixReal& m) const;
        
    };
    
}



#endif /* NeighborJoining_hpp */
