#ifndef TreeBipartitions_H
#define TreeBipartitions_H

#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <boost/dynamic_bitset.hpp>


#include <vector>
#include <string>

namespace RevBayesCore {
    
    class TreeBipartitions : public TypedFunction< RbVector<boost::dynamic_bitset<> > > {
        
    public:
        TreeBipartitions(const TypedDagNode<Tree> *t);                                                  //!< Default constructor
        virtual      ~TreeBipartitions(void);                                                           //!< Destructor
        
        // Basic utility functions
        TreeBipartitions*                           clone(void) const;                                  //!< Clone object
        void                                        update(void);                                       //!< Clone the function
        const std::vector<double>&                  getBipartitionAges(void);
        
    protected:
        void swapParameterInternal(const DagNode *oldP, const DagNode *newP);                           //!< Implementation of swaping parameters
        
    private:
        // members
        const TypedDagNode<Tree>*                   tree;
        std::vector<boost::dynamic_bitset<> >       bipartitions;
        std::map<boost::dynamic_bitset<>, const TopologyNode*> bipartitionMap;
        std::vector<double> bipartitionAges;
        
        void computeBipartitions();
        void computeBipartitions(const TopologyNode* node, std::map <const TopologyNode*, unsigned long>& nodeToBitVectorIndex, const std::map<std::string, size_t>& nameToIndex );

    };
    
}

#endif

