#ifndef UPGMA_hpp
#define UPGMA_hpp

#include <cstdio>
#include <vector>
#include <set>

namespace RevBayesCore {
class Clade;
class DistanceMatrix;
class Taxon;
class TopologyNode;
class Tree;
    
    /**
     * This class provides the UPGMA algorithm to construct a tree from a distance matrix.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2013-04-15, version 1.0
     */
    class UPGMA {
        
    public:
        
//        enum SIM_CONDITION { TIME, SURVIVAL, ROOT };
        
        UPGMA();
        
        Tree*                   constructTree( const DistanceMatrix& d ) const;
        
    private:
        
        TopologyNode*           constructTreeRecursively( std::vector<TopologyNode*> &nodes, MatrixReal& m) const;
        
    };
    
}



#endif /* UPGMA_hpp */
