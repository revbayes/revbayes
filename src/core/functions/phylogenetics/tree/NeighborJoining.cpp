#include "DistanceMatrix.h"
#include "RbConstants.h"
#include "Taxon.h"
#include "TopologyNode.h"
#include "Tree.h"
#include "NeighborJoining.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"


using namespace RevBayesCore;


NeighborJoining::NeighborJoining( void )
{
    
}


Tree* NeighborJoining::constructTree(const DistanceMatrix &d) const
{
    // get some information about the data
    const std::vector<Taxon>&   taxa        = d.getTaxa();
    size_t                      num_tips    = d.getSize();
    MatrixReal                  distances   = d.getMatrix();
    
    // make sure the diagonal is set to 0.0
    for ( size_t i=0; i<num_tips; ++i )
    {
        distances[i][i] = 0.0;
    }
    
    // first, we need to create a vector with all the nodes
    std::vector<TopologyNode*> active_nodes = std::vector<TopologyNode*>(num_tips, NULL);
    for (size_t i=0; i<num_tips; ++i)
    {
        TopologyNode* node = new TopologyNode( taxa[i], i );
        active_nodes[i] = node;
    }
    
    // perform the NeighborJoining algorithm recursively
    TopologyNode* root = constructTreeRecursively(active_nodes, distances);
    
    // construct the tree
    Tree* nj_tree = new Tree();
    nj_tree->setRoot(root, true);
        
    // finally, return our constructed tree
    return nj_tree;
}



TopologyNode* NeighborJoining::constructTreeRecursively(std::vector<TopologyNode *> &active_nodes, MatrixReal &distances) const
{
    
    size_t n = distances.size();
    
    // compute the Q matrix
    std::vector<double> row_sums = std::vector<double>(n, 0.0);
    for (size_t i=0; i<n; ++i)
    {
        for (size_t j=0; j<n; ++j)
        {
            if ( i != j )
            {
                row_sums[i] += distances[i][j];
            }
        }
    }
    
    double min_Q = RbConstants::Double::inf;
    std::set< std::pair<size_t,size_t> > min_pairs;
//    MatrixReal Q = MatrixReal(n);
    for (size_t i=0; i<n; ++i)
    {
        
        for (size_t j=i+1; j<n; ++j)
        {
            double this_Q = (n-2)*distances[i][j] - row_sums[i] - row_sums[j];
//            Q[i][j] = this_Q;
            if ( this_Q < min_Q )
            {
                // set the new minimum
                min_Q = this_Q;
                
                // clear the previously best entries
                min_pairs.clear();
                
                // enter the new entry i <=> j
                min_pairs.insert( std::pair<size_t,size_t>(i,j) );

            }
            else if ( this_Q == min_Q )
            {
                
                // enter the new entry i <=> j
                min_pairs.insert( std::pair<size_t,size_t>(i,j) );
            }
            
        }
        
    }
    
    // find the smallest distances
    std::pair<size_t,size_t> best_pair;
    size_t num_best_pairs = min_pairs.size();
    if ( num_best_pairs == 1 )
    {
        best_pair = *min_pairs.begin();
    }
    else
    {
        // Get random number generator
        RandomNumberGenerator* rng = GLOBAL_RNG;
        double u = rng->uniform01();
        std::set< std::pair<size_t,size_t> >::iterator it = min_pairs.begin();
        while ( u > (1.0/num_best_pairs) )
        {
            u -= (1.0/num_best_pairs);
            ++it;
        }
        best_pair = *it;
       
    }
    size_t index_A = best_pair.first;
    size_t index_B = best_pair.second;
    
    
    // make sure that A < B
    if ( index_B < index_A )
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
    
    // set the branch lengths
    double bl_left = distances[index_A][index_B] / 2.0 + 1.0/(2.0*(n-2.0)) * ( row_sums[index_A] - row_sums[index_B] );
    double bl_right = distances[index_A][index_B] - bl_left;
    left->setBranchLength( bl_left );
    right->setBranchLength( bl_right );

    // erase the two children from the vector of active nodes
    active_nodes.erase( active_nodes.begin()+index_B );
    active_nodes.erase( active_nodes.begin()+index_A );
    
    // add the new parent at the end
    active_nodes.push_back( parent );
    
    // create a new row and column in the distance matrix
    size_t num_elements = distances.size();
    
    distances.addColumn();
    distances.addRow();
    // ensure that the diagonal stays infinite
    distances[num_elements][num_elements] = RbConstants::Double::inf;
    for (size_t i=0; i<num_elements; ++i)
    {
        if ( i != index_A && i != index_B )
        {
            double new_dist = (distances[index_A][i] + distances[index_B][i] - distances[index_A][index_B]) / 2.0;
            distances[num_elements][i] = new_dist;
            distances[i][num_elements] = new_dist;
        }
    }
    
    // delete the old rows and columns for A and B
    distances.deleteRow(index_B);
    distances.deleteRow(index_A);
    distances.deleteColumn(index_B);
    distances.deleteColumn(index_A);
    
    
    // finally, continue the recursive call
    if ( num_elements > 4 )
    {
        return constructTreeRecursively(active_nodes, distances);
    }
    else
    {
        
        // get the corresponding nodes
        TopologyNode* first_node  = active_nodes[0];
        TopologyNode* second_node = active_nodes[1];
        TopologyNode* third_node  = active_nodes[2];
        
        // create the new parent node
        TopologyNode* root = new TopologyNode();
        
        // join the two nodes
        root->addChild( first_node );
        root->addChild( second_node );
        root->addChild( third_node );
        first_node->setParent( root );
        second_node->setParent( root );
        third_node->setParent( root );
        
//        double bl_first = distances[0][1]/2.0 + (distances[0][1] + distances[0][2] - distances[1][0] - distances[1][2] )/2.0
        double bl_first  = distances[0][1]/2.0 + (distances[0][2] - distances[1][2] )/2.0;
        double bl_second = distances[0][1] - bl_first;
        double bl_third  = distances[0][2] - bl_first;
        
        first_node->setBranchLength( bl_first );
        second_node->setBranchLength( bl_second );
        third_node->setBranchLength( bl_third );

        return root;
    }
}
