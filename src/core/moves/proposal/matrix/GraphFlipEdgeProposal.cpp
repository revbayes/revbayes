#include "GraphFlipEdgeProposal.h"

#include <iostream>
#include <vector>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbVector.h"
#include "Cloneable.h"
#include "MatrixReal.h"
#include "StochasticNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
GraphFlipEdgeProposal::GraphFlipEdgeProposal( StochasticNode<MatrixReal> *n, const RbVector<RbVector<std::int64_t> >& e, double l, bool s) : Proposal(),
array(NULL),
matrix( n ),
edges( e ),
sampling_probability( l ),
storedValue( 0.0 ),
symmetric( s )
{
    // build the complete edge list if no edges are provided
    edge_list_length = edges.size();
    if (edge_list_length == 0) {
        
        for (size_t i = 0; i < n->getValue().size(); i++) {
            for (size_t j = i+1; j < n->getValue().size(); j++) {
                std::vector<std::int64_t> tmp;
                tmp.push_back(i+1);
                tmp.push_back(j+1);
                edges.push_back( tmp );
            }
        }
        
        edge_list_length = edges.size();
    }
    
    // tell the base class to add the node
    addNode( matrix );
    
}

/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void GraphFlipEdgeProposal::cleanProposal( void )
{
    // clear touched edge elements set
    touched_edge_elements = std::set<size_t>();
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
GraphFlipEdgeProposal* GraphFlipEdgeProposal::clone( void ) const
{
    
    return new GraphFlipEdgeProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& GraphFlipEdgeProposal::getProposalName( void ) const
{
    static std::string name = "GraphFlipEdgeMove";
    
    return name;
}


double GraphFlipEdgeProposal::getProposalTuningParameter( void ) const
{
    return sampling_probability;
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
double GraphFlipEdgeProposal::doProposal( void )
{
    
    // Get random number generator
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    // clear touched edge elements set
    touched_edge_elements = std::set<size_t>();
    
    // sample the first edge to update w.p. 1
    size_t first_index = size_t( rng->uniform01() * edge_list_length );
    touched_edge_elements.insert( first_index );
    
    // sample additional edges to update w.p. p
    if (sampling_probability > 0.0) {
        for (size_t i = 0; i < edge_list_length; i++) {
            if (rng->uniform01() < sampling_probability) {
                touched_edge_elements.insert( i );
            }
        }
    }
    
    // flip and touch all sampled edges
    MatrixReal& v = matrix->getValue();
    for (std::set<size_t>::iterator it = touched_edge_elements.begin(); it != touched_edge_elements.end(); it++) {
        size_t i = edges[ *it ][ 0 ] - 1;
        size_t j = edges[ *it ][ 1 ] - 1;
        v[i][j] = ( v[i][j] == 1 ? 0 : 1 );
        if (symmetric) {
            v[j][i] = v[i][j];
            matrix->addTouchedElementIndex( j * v.getNumberOfRows() + i );
        }
        matrix->addTouchedElementIndex( i * v.getNumberOfRows() + j );

    }
    
//    if (!v.isSymmetric()) {
//        throw RbException("matrix is not symmetric!\n");
//    }

    // symmetric proposal
    return 0;
}


/**
 *
 */
void GraphFlipEdgeProposal::prepareProposal( void )
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
void GraphFlipEdgeProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
    o << "p = ";
    if (name_only == false)
    {
        o << sampling_probability;
    }
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void GraphFlipEdgeProposal::undoProposal( void )
{

    // revert all flipped edges
    MatrixReal& v = matrix->getValue();
    for (std::set<size_t>::iterator it = touched_edge_elements.begin(); it != touched_edge_elements.end(); it++) {
        size_t i = edges[ *it ][ 0 ] - 1;
        size_t j = edges[ *it ][ 1 ] - 1;
        v[i][j] = ( v[i][j] == 1 ? 0 : 1 );
        if (symmetric) {
            v[j][i] = v[i][j];
        }
    }
    
    // clear all touched elements
    matrix->clearTouchedElementIndices();
    
//    if (!v.isSymmetric()) {
//        throw RbException("matrix is not symmetric!\n");
//    }
    
    // clear touched edge elements set
    touched_edge_elements = std::set<size_t>();
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void GraphFlipEdgeProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    if (oldN == array)
    {
        array = static_cast< StochasticNode<RbVector<RbVector<double> > >* >(newN) ;
    }
    else if (oldN == matrix)
    {
        matrix = static_cast< StochasticNode<MatrixReal>* >(newN) ;
    }
    
}


void GraphFlipEdgeProposal::setProposalTuningParameter(double tp)
{
    sampling_probability = tp;
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
void GraphFlipEdgeProposal::tune( double rate )
{
    
    if ( rate > 0.44 )
    {
        sampling_probability *= (1.0 + ((rate-0.44)/0.56) );
    }
    else
    {
        sampling_probability /= (2.0 - rate/0.44 );
    }
    
    // bound value s.t. 0 < switch_probability <= 1.0
    if (sampling_probability > 1.0)
    {
        sampling_probability = 1.0;
    }
    else if (sampling_probability < 1.0 / edge_list_length)
    {
        sampling_probability = 1.0 / edge_list_length;
    }
}

