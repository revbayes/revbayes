#include "DistributionUniform.h"
#include "NodeTimeSlideUniformCharacterHistoryProposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbException.h"
#include "TreeUtilities.h"
#include "TypedDagNode.h"

#include <cmath>
#include <iostream>

using namespace RevBayesCore;

///**
// * Constructor
// *
// * Here we simply allocate and initialize the Proposal object.
// */
//NodeTimeSlideUniformCharacterHistoryProposal::NodeTimeSlideUniformCharacterHistoryProposal( StochasticNode<Tree> *n ) : Proposal(),
//    variable( n )
//{
//    // tell the base class to add the node
//    addNode( variable );
//    
//}
//
//
///**
// * The cleanProposal function may be called to clean up memory allocations after AbstractMove
// * decides whether to accept, reject, etc. the proposed value.
// *
// */
//void NodeTimeSlideUniformCharacterHistoryProposal::cleanProposal( void )
//{
//    ; // do nothing
//}
//
///**
// * The clone function is a convenience function to create proper copies of inherited objected.
// * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
// *
// * \return A new copy of the proposal.
// */
//NodeTimeSlideUniformCharacterHistoryProposal* NodeTimeSlideUniformCharacterHistoryProposal::clone( void ) const
//{
//    
//    return new NodeTimeSlideUniformCharacterHistoryProposal( *this );
//}
//
//
///**
// * Get Proposals' name of object
// *
// * \return The Proposals' name.
// */
//const std::string& NodeTimeSlideUniformCharacterHistoryProposal::getProposalName( void ) const
//{
//    static std::string name = "NodeTimeSlideUniformCharacterHistory";
//    
//    return name;
//}
//
//
///**
// * Perform the proposal.
// *
// * A Uniform-simplex proposal randomly changes some values of a simplex, although the other values
// * change too because of the renormalization.
// * First, some random indices are drawn. Then, the proposal draws a new somplex
// *   u ~ Uniform(val[index] * alpha)
// * where alpha is the tuning parameter.The new value is set to u.
// * The simplex is then renormalized.
// *
// * \return The hastings ratio.
// */
//double NodeTimeSlideUniformCharacterHistoryProposal::doProposal( void )
//{
//    
//    // Get random number generator
//    RandomNumberGenerator* rng     = GLOBAL_RNG;
//    
//    Tree& tau = variable->getValue();
//    
//    // pick a random node which is not the root and neithor the direct descendant of the root
//    TopologyNode* node;
//    do {
//        double u = rng->uniform01();
//        size_t index = size_t( std::floor(tau.getNumberOfNodes() * u) );
//        node = &tau.getNode(index);
//    } while ( node->isRoot() || node->isTip() );
//    
//    TopologyNode& parent = node->getParent();
//    
//    // we need to work with the times
//    double parent_age  = parent.getAge();
//    double my_age      = node->getAge();
//    double child_Age   = node->getChild( 0 ).getAge();
//    if ( child_Age < node->getChild( 1 ).getAge())
//    {
//        child_Age = node->getChild( 1 ).getAge();
//    }
//    
//    // now we store all necessary values
//    storedNode = node;
//    storedAge = my_age;
//    
//    // draw new ages and compute the hastings ratio at the same time
//    double my_new_age = (parent_age-child_Age) * rng->uniform01() + child_Age;
//    
//    // set the age
//    tau.getNode(node->getIndex()).setAge( my_new_age );
//    
//    return 0.0;
//    
//}
//
//
///**
// *
// */
//void NodeTimeSlideUniformCharacterHistoryProposal::prepareProposal( void )
//{
//    
//}
//
//
///**
// * Print the summary of the Proposal.
// *
// * The summary just contains the current value of the tuning parameter.
// * It is printed to the stream that it passed in.
// *
// * \param[in]     o     The stream to which we print the summary.
// */
//void NodeTimeSlideUniformCharacterHistoryProposal::printParameterSummary(std::ostream &o) const
//{
//    
//}
//
//
///**
// * Reject the Proposal.
// *
// * Since the Proposal stores the previous value and it is the only place
// * where complex undo operations are known/implement, we need to revert
// * the value of the variable/DAG-node to its original value.
// */
//void NodeTimeSlideUniformCharacterHistoryProposal::undoProposal( void )
//{
//    
//    // undo the proposal
//    storedNode->setAge( storedAge );
//    
//}
//
//
///**
// * Swap the current variable for a new one.
// *
// * \param[in]     oldN     The old variable that needs to be replaced.
// * \param[in]     newN     The new RevVariable.
// */
//void NodeTimeSlideUniformCharacterHistoryProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
//{
//    
//    variable = static_cast<StochasticNode<Tree>* >(newN) ;
//    
//}
//
//
///**
// * Tune the Proposal to accept the desired acceptance ratio.
// *
// * The acceptance ratio for this Proposal should be around 0.44.
// * If it is too large, then we increase the proposal size,
// * and if it is too small, then we decrease the proposal size.
// */
//void NodeTimeSlideUniformCharacterHistoryProposal::tune( double rate )
//{
//    
//}

