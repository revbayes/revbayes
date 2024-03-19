#include <cstddef>
#include <vector>

#include "DagNodeMap.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;



DagNode*& DagNodeMap::operator[]( const DagNode *d )
{
    
    size_t s=sources.size();
    for (size_t i=0; i<s; ++i)
    {
        
        if ( sources[i] == d )
        {
            return targets[i];
        }
        
    }
    
    // the node wasn't there
    sources.push_back( d );
    targets.push_back( NULL );
    return targets[targets.size()-1];
}


std::vector<DagNode*>::iterator DagNodeMap::begin()
{
    return targets.begin();
}


std::vector<DagNode*>::iterator DagNodeMap::find(const DagNode* n)
{
    size_t s=sources.size();
    for (size_t i=0; i<s; ++i)
    {
        
        if ( sources[i] == n )
        {
            return targets.begin() + i;
        }
        
    }
    
    return targets.end();
}


std::vector<DagNode*>::iterator DagNodeMap::end()
{
    return targets.end();
}




