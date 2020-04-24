#include "DistanceMatrix.h"
#include "Taxon.h"
#include "TopologyNode.h"
#include "Tree.h"
#include "UPGMA.h"


using namespace RevBayesCore;


UPGMA::UPGMA( void )
{
    
}


Tree* UPGMA::constructTree(const DistanceMatrix &d) const
{
    // get some information about the data
    const std::vector<Taxon>&   taxa        = d.getTaxa();
    size_t                      num_tips    = d.getSize();
    MatrixReal                  distances   = d.getMatrix();
    
    // first, we need to create a vector with all the nodes
    std::vector<TopologyNode*> active_nodes = std::vector<TopologyNode*>(num_tips, NULL);
    for (size_t i=0; i<num_tips; ++i)
    {
        TopologyNode* node = new TopologyNode( taxa[i], i );
        node->setAge( taxa[i].getAge() );
        active_nodes[i] = node;
    }
    
    // perform the UPGMA algorithm recursively
    TopologyNode* root = constructTreeRecursively(active_nodes, distances);
    
    // construct the tree
    Tree* upgma_tree = new Tree();
    upgma_tree->setRoot(root, false);
    
    // finally, return our constructed tree
    return upgma_tree;
}



TopologyNode* UPGMA::constructTreeRecursively(std::vector<TopologyNode *> &active_nodes, MatrixReal &distances) const
{
    
    // find the smallest distances
    size_t index_A = 0;
    size_t index_B = 1;
    distances.getIndexOfMin(index_A, index_B);
    
    // make sure that A < B
    if ( index_B > index_A )
    {
        size_t tmp = index_A;
        index_A = index_B;
        index_B = tmp;
    }
    
    // get the corresponding nodes
    TopologyNode* left  = active_nodes[index_A];
    TopologyNode* right = active_nodes[index_B];
    
    // create the new parent node
    TopologyNode* parent = new TopologyNode();
    
    // join the two nodes
    parent->addChild( left );
    parent->addChild( right );
    left->setParent( parent );
    right->setParent( parent );
    
    // set the age of the parent
    double parent_age = distances[index_A][index_B] / 2.0;
    parent->setAge( parent_age );
    
    // erase the two children from the vector of active nodes
    active_nodes.erase( active_nodes.begin()+index_B );
    active_nodes.erase( active_nodes.begin()+index_A );
    
    // add the new parent at the end
    active_nodes.push_back( parent );
    
    // create a new row and column in the distance matrix
    size_t num_elements = distances.size();
    
    distances.addColumn();
    distances.addRow();
    for (size_t i=0; i<num_elements; ++i)
    {
        if ( i != index_A && i != index_B )
        {
            double new_dist = (distances[index_A][i] + distances[index_B][i]) / 2.0;
            distances[num_elements][i] = new_dist;
            distances[i][num_elements] = new_dist;
        }
    }
    
    // delete the old rows and columns for A and B
    distances.deleteRow(index_A);
    distances.deleteRow(index_B);
    distances.deleteColumn(index_A);
    distances.deleteColumn(index_B);
    
    
    // finally, continue the recursive call
    if ( num_elements > 2 )
    {
        return constructTreeRecursively(active_nodes, distances);
    }
    else
    {
        return parent;
    }
}
