#include "ShortestDistanceFunction.h"
#include <set>
#include <vector>
#include <cstddef>
#include <cstdint>
#include "RbConstants.h"
#include "RbVector.h"
#include "TypedFunction.h"
#include "Cloner.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;


/**
 * Default Constructor for the ShortestDistanceFunction
 *
 * @param adj a vector of vectors of longs that represent a matrix denoting the connectivity of nodes. adj[i][j] > 0 for connected nodes and 0 for disconnected nodes.
 * @param dist a vector of vectors of doubles that represent a matrix denoting the distance between connected nodes. dist[i][j] denotes the distance between nodes i and j.
 */

ShortestDistanceFunction::ShortestDistanceFunction(const TypedDagNode<RbVector<RbVector<std::int64_t> > >* adj, const TypedDagNode<RbVector<RbVector<double> > >* dist) : TypedFunction<RbVector<RbVector<double> > >( new RbVector<RbVector<double> >() ),
num_nodes( adj->getValue().size() ),
adjacencies( adj ),
distances( dist )

{
    
    *value = RbVector<RbVector<double> >(num_nodes, RbVector<double>(num_nodes, 1));
   
    addParameter( adjacencies );
    addParameter( distances );
    
    update();
}


ShortestDistanceFunction* ShortestDistanceFunction::clone( void ) const
{
    return new ShortestDistanceFunction(*this);
}

/**
 * ShortestDistanceFunction::CreateAdjacencySets
 *
 * Creates a set for each node that denotes which nodes it is connected to
 *
 * @param adj a vector of vectors of type long
 *
 * @return s a vector of sets where each set denotes the adjacent nodes for that index in the vector
 *
 */
std::vector<std::set<size_t> > ShortestDistanceFunction::createAdjacencySets(const RbVector<RbVector<std::int64_t> >& adj)
{

    std::vector<std::set<size_t> > s;
    for (size_t i = 0; i < num_nodes; i++)
    {
        std::set<size_t> t;
        for (size_t j = 0; j < num_nodes; j++)
        {
            if (i != j && adj[i][j] > 0)
            {
                t.insert(j);
            }
        }
        s.push_back(t);
    }
    
    return s;
}

/**
 * ShortestDistanceFunction::findShortestPaths
 *
 * Runs Dijkstra's shortest path algorithm
 *
 * @param adj
 * @param dist
 *
 */

RbVector<RbVector<double> > ShortestDistanceFunction::findShortestPaths(const RbVector<RbVector<std::int64_t> >& adj, const RbVector<RbVector<double> >& dist)
{
    
    std::vector<std::set<size_t> > adj_set = createAdjacencySets(adj);
    
    // 0. Keep track of final distances for all origin nodes under Dijkstra's
    RbVector<RbVector<double> > return_distances(num_nodes, std::vector<double>(num_nodes, RbConstants::Double::inf));
    
    // run Dijkstra's for each origin node
    for (size_t origin_node = 0; origin_node < num_nodes; origin_node++)
    {
        // 1. Set the initial node as current.
        size_t current_node = origin_node;
        
        // 2. Assign to every node a tentative distance value: set it to zero for our initial node and to infinity for all other nodes.
        return_distances[current_node][current_node] = 0.0;

        // Mark all other nodes unvisited.
        std::vector<bool> visited_list(num_nodes, false);
        visited_list[current_node] = true;
        
        // Create a set of all the unvisited nodes called the unvisited set.
        std::set<size_t> unvisited_nodes;
        for (size_t j = 0; j < num_nodes; j++)
        {
            if (j != current_node) {
                unvisited_nodes.insert(j);
            }
        }

        bool stop = false;
        while (!stop)
        {
        
            // 3. For the current node, consider all of its unvisited neighbors and calculate their tentative distances.
            //    Compare the newly calculated tentative distance to the current assigned value and assign the smaller one.
            //    For example, if the current node A is marked with a distance of 6, and the edge connecting it with a neighbor B has length 2, then the distance to B (through A) will be 6 + 2 = 8.
            //    If B was previously marked with a distance greater than 8 then change it to 8. Otherwise, keep the current value.
        
//            std::cout << "\t" << current_node << "\n";
            for (std::set<size_t>::iterator it = adj_set[current_node].begin(); it != adj_set[current_node].end(); it++)
            {
                size_t j = *it;
                
                // do not revisit visited nodes
                if (visited_list[j] == true)
                    continue;
                
//                std::cout << "\t\t" << j << "\n";
//                std::cout << "\t\tcurr_dist " << return_distances[origin_node][j] << "\n";
//                std::cout << "\t\tposs_dist " << return_distances[origin_node][current_node] + dist[current_node][j] << "\n";
                if (return_distances[origin_node][current_node] + dist[current_node][j] < return_distances[origin_node][j])
                {
                    return_distances[origin_node][j] = return_distances[origin_node][current_node] + dist[current_node][j];
                }
//                 std::cout << "\t\tfinal_dist " << return_distances[origin_node][j] << "\n";
            }
            
            
            // 4. When we are done considering all of the neighbors of the current node, mark the current node as visited and remove it from the unvisited set. A visited node will never be checked again.
            visited_list[current_node] = true;
            unvisited_nodes.erase(current_node);
            
            
            // 5. If the smallest tentative distance among the nodes in the unvisited set is infinity (when planning a complete traversal; occurs when there is no connection between the initial node and remaining unvisited nodes), then stop. The algorithm has finished.
            double smallest_value = RbConstants::Double::inf;
            bool found_smallest = false;
            for (std::set<size_t>::iterator jt = unvisited_nodes.begin(); jt != unvisited_nodes.end(); jt++)
            {

                size_t j = *jt;

                // the next current_node has the smallest distance
                if (return_distances[origin_node][j] < smallest_value)
                {
                    found_smallest = true;
                    smallest_value = return_distances[origin_node][j];
                    current_node = j;
                }
            }
        
            // all connected nodes have been evaluated
            if (!found_smallest)
            {
                stop = true;
            }
            
        }
    }
    
    return return_distances;
    
}

void ShortestDistanceFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == adjacencies)
    {
        adjacencies = static_cast<const TypedDagNode<RbVector<RbVector<std::int64_t> > >* >( newP );
    }
    else if (oldP == distances)
    {
        distances = static_cast<const TypedDagNode<RbVector<RbVector<double> > >* >( newP );
    }
    
}

void ShortestDistanceFunction::update( void )
{
    // get the new values
    const RbVector<RbVector<std::int64_t> >& adj = adjacencies->getValue();
    const RbVector<RbVector<double> >& dist = distances->getValue();
    
    // find all the shortest paths
    *value = findShortestPaths(adj, dist);
    
}
