#include <cmath>
#include <iostream>
#include <cstddef>
#include <set>
#include <type_traits>
#include <vector>

#include "DistributionBinomial.h"
#include "DistributionPoisson.h"
#include "GraphFlipCliqueProposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbVector.h"
#include "MatrixReal.h"
#include "Proposal.h"
#include "RbVectorImpl.h"
#include "StochasticNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
GraphFlipCliqueProposal::GraphFlipCliqueProposal( StochasticNode<MatrixReal> *n, const RbVector<std::int64_t>& v, double sr, double vp, double ep, bool s) : Proposal(),
array(NULL),
matrix( n ),
vertices( v ),
set_sample_rate( sr ),
vertex_probability( vp ),
edge_probability( ep ),
storedValue( 0.0 ),
symmetric( s ),
undo_needed( false )
{
    // build the complete edge list if no edges are provided
    vertex_list_length = vertices.size();
    if (vertex_list_length == 0) {
        
        for (size_t i = 0; i < n->getValue().size(); i++) {
            vertices.push_back(i+1);
        }
        vertex_list_length = vertices.size();
    }
    
    if (edge_probability > 0.5) edge_probability = 0.5;
    
    // tell the base class to add the node
    addNode( matrix );
    
}

/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void GraphFlipCliqueProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
GraphFlipCliqueProposal* GraphFlipCliqueProposal::clone( void ) const
{
    
    return new GraphFlipCliqueProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& GraphFlipCliqueProposal::getProposalName( void ) const
{
    static std::string name = "GraphFlipCliqueMove";
    
    return name;
}


double GraphFlipCliqueProposal::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}


/**
 * Perform the proposal.
 *
 * A sliding proposal draws a random uniform number u ~ unif (-0.5,0.5)
 * and MatrixRealSingleElementSlidings the current vale by
 * delta = lambda * u
 * where lambda is the tuning parameter of the proposal to influence the size of the proposals.
 *
 * \return The hastings ratio.
 */
double GraphFlipCliqueProposal::doProposal( void )
{
    
    // to be returned
    double lnProposalRatio = 0.0;
    double ln_p = log(edge_probability);
    double ln_q = log(1.0 - edge_probability);
    
    // undo not necessarily needed
//    undo_needed = false;
    
    // clear touched edge vector
    stored_edge_elements = std::vector<Edge>();
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    MatrixReal& v = matrix->getValue();
    
    std::cout << "before\n";
    for (size_t a = 0; a < v.size(); a++) {
        for (size_t b = 0; b < v.size(); b++) {
            std::cout << ( v[a][b] == 0 ? "  " : ". " );
        }
        std::cout << "\n";
    }
    std::cout << "\n";
    
    
    // Always flip at least one clique
    size_t n_to_shift = 1 + RbStatistics::Poisson::rv( set_sample_rate-1.0, *GLOBAL_RNG );
    
    std::cout << "n_to_shift = " << n_to_shift << "\n";
    
    // repeat the FlipClique procedure multiple times
    for (size_t i = 0; i < n_to_shift; i++) {
        
        // Are we adding or deleting edges within a vertex-set?
        size_t new_edge_state = ( rng->uniform01() < 0.5 ? 0 : 1 );
        size_t other_edge_state = ( new_edge_state == 0 ? 1 : 0 );
        
        // Sample set of nodes in potential clique
        std::set<size_t> unique_possible_clique_vertices;
        
        // clique must contain at least two nodes
        size_t max_clique_size = 2 + RbStatistics::Binomial::rv( vertices.size()-2, vertex_probability, *GLOBAL_RNG );
        std::cout << "sampled vertex for set " << i << " = { ";
        for (size_t j = 0; j < max_clique_size; j++)
        {
            size_t idx = size_t(vertices.size() * rng->uniform01());
            std::cout << idx << " ";
            unique_possible_clique_vertices.insert( idx );
        }
        std::cout << "}\n";
        
        // convert set to vertex for iteration
        std::vector<size_t> possible_clique_vertices( unique_possible_clique_vertices.begin(), unique_possible_clique_vertices.end() );
        
        // sample edges w.p. edge_probability between each pair of vertices
        std::cout << "sampled edges = {\n";
        std::vector<size_t> num_constant(2, 0);
//        std::vector<size_t> new_num_edges(2, 0);
        for (size_t i = 0; i < possible_clique_vertices.size(); i++)
        {
            size_t idx = vertices[ possible_clique_vertices[i] ]-1;
            for (size_t j = i+1; j < possible_clique_vertices.size(); j++)
            {
                size_t jdx = vertices[ possible_clique_vertices[j] ]-1;
                
                // for proposal ratio
//                old_num_edges[ v[idx][jdx] ] += 1;
                
                // do we flip the edge to new_edge_state?
                double u = rng->uniform01();
                if (v[idx][jdx] == new_edge_state)
                {
                    std::cout << "\t* " << idx << "," << jdx << " : " << v[idx][jdx] << " -> " << new_edge_state << "\n";
                    num_constant[new_edge_state] += 1;
                }
                else if (u < edge_probability)
                {
                    stored_edge_elements.push_back( Edge(idx, jdx, v[idx][jdx]) );
                    std::cout << "\t+ " << idx << "," << jdx << " : " << v[idx][jdx] << " -> " << new_edge_state << "\n";
                    
                    v[idx][jdx] = new_edge_state;
                    matrix->addTouchedElementIndex( idx * v.getNumberOfRows() + jdx );
                    if (symmetric) {
                        v[jdx][idx] = v[idx][jdx];
                        matrix->addTouchedElementIndex( jdx * v.getNumberOfRows() + idx );
                    }
                    
                }
                else
                {
                    std::cout << "\t- " << idx << "," << jdx << " : " << v[idx][jdx] << " -> " << other_edge_state << "\n";
                    if (v[idx][jdx] == other_edge_state) {
                        num_constant[other_edge_state] += 1;
                    }
                }
                
                // for proposal ratio
//                new_num_edges[ v[idx][jdx] ] += 1;
                
            }
        }
        std::cout << "}\n";
        
        // Example:
        // 00011 -> 01111; Prob( (3,2)->(1,4) | p ) = Prob(n->n' | p) = q * p * p * 1 * 1 = p^2 * q
        // 01111 -> 00011; Prob( (1,4)->(3,2) | p ) = Prob(n'->n | p) = 1 * p * p * q * q = p^2 * q^2
        
        // update proposal ratio
//        lnProposalRatio += ln_p * old_num_edges[ other_edge_state ] + ln_q * old_num_edges[ new_edge_state ]; // reverse move
//        lnProposalRatio -= ln_p * new_num_edges[ new_edge_state ] + ln_q * new_num_edges[ other_edge_state ]; // forward move
        
        // p( x' -> x ) =
        // p( x -> x' ) =
//        lnProposalRatio -= ln_q * old_num_edges[ new_edge_state ];
        
        std::cout << "new_edge_state   = " << new_edge_state << "\n";
        std::cout << "other_edge_state = " << other_edge_state << "\n";
        std::cout << "n = " << num_constant[0] << "," << num_constant[1] << "\n";
        lnProposalRatio += ln_q * (num_constant[other_edge_state] - num_constant[new_edge_state]);
    }
    
    std::cout << "after\n";
    for (size_t a = 0; a < v.size(); a++) {
        for (size_t b = 0; b < v.size(); b++) {
            std::cout << ( v[a][b] == 0 ? "  " : ". " );
        }
        std::cout << "\n";
    }
    std::cout << "\n";
    
    return lnProposalRatio;
}


/**
 *
 */
void GraphFlipCliqueProposal::prepareProposal( void )
{
    
}


/**
 * Print the summary of the Proposal.
 *
 * The summary just contains the current value of the tuning parameter.
 * It is printed to the stream that it passed in.
 *
 * \param[in]     o     The stream to which we print the summary.
 */
void GraphFlipCliqueProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
    o << "set_sample_rate = " << set_sample_rate;
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void GraphFlipCliqueProposal::undoProposal( void )
{
    
    MatrixReal& v = matrix->getValue();
    
    std::cout << "pre-undo\n";
    for (size_t a = 0; a < v.size(); a++) {
        for (size_t b = 0; b < v.size(); b++) {
            std::cout << ( v[a][b] == 0 ? "  " : ". " );
        }
        std::cout << "\n";
    }
    std::cout << "\n";
    
    for (std::vector<Edge>::reverse_iterator rit = stored_edge_elements.rbegin(); rit != stored_edge_elements.rend(); rit++)
    {
        v[ rit->from ][ rit->to ] = rit->value;
        //        std::cout << "UNDO " << rit->from << "," << rit->to << " = " << rit->value << "\n";
        if (symmetric) {
            v[ rit->to ][ rit->from ] = rit->value;
        }
    }
    
    // clear all touched elements
    matrix->clearTouchedElementIndices();
    
    std::cout << "post-undo\n";
    for (size_t a = 0; a < v.size(); a++) {
        for (size_t b = 0; b < v.size(); b++) {
            std::cout << ( v[a][b] == 0 ? "  " : ". " );
        }
        std::cout << "\n";
    }
    std::cout << "\n";
    
    
    // clear touched edge elements set
    stored_edge_elements = std::vector<Edge>();
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void GraphFlipCliqueProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    if (oldN == array)
    {
        array = static_cast< StochasticNode<RbVector<RbVector<double> > >* >(newN) ;
    }
    else
    {
        matrix = static_cast< StochasticNode<MatrixReal>* >(newN) ;
    }
    
}


void GraphFlipCliqueProposal::setProposalTuningParameter(double tp)
{
    // this proposal has no tuning parameter: nothing to do
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
void GraphFlipCliqueProposal::tune( double rate )
{
//    
//    if ( rate > 0.44 )
//    {
//        sampling_probability *= (1.0 + ((rate-0.44)/0.56) );
//    }
//    else
//    {
//        sampling_probability /= (2.0 - rate/0.44 );
//    }
//    
//    // bound value s.t. 0 < switch_probability <= 1.0
//    if (sampling_probability > 1.0)
//    {
//        sampling_probability = 1.0;
//    }
//    else if (sampling_probability < 1.0 / vertex_list_length)
//    {
//        sampling_probability = 1.0 / vertex_list_length;
//    }
}

