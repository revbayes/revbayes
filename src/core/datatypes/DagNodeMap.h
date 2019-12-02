#ifndef DagNodeMap_H
#define DagNodeMap_H

#include <vector>

namespace RevBayesCore {
    
    class DagNode;
    
    class DagNodeMap {
        
    public:
        virtual                                 ~DagNodeMap(void) {}
        
        DagNode*&                               operator[]( const DagNode *d );

        std::vector<DagNode*>::iterator         begin();
        std::vector<DagNode*>::iterator         find(const DagNode* n);
        std::vector<DagNode*>::iterator         end();
        
    private:
        
        std::vector<const DagNode*>             sources;
        std::vector<DagNode*>                   targets;
        
    };
    
}

#endif

