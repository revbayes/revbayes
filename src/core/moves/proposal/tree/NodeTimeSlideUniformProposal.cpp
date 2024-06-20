#include <cstddef>
#include <cmath>
#include <iostream>

#include "DistributionUniform.h"
#include "NodeTimeSlideUniformProposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "Proposal.h"
#include "StochasticNode.h"
#include "TopologyNode.h"
#include "Tree.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
NodeTimeSlideUniformProposal::NodeTimeSlideUniformProposal( StochasticNode<Tree> *n, StochasticNode< RbVector<Tree> >* vec_n ) : Proposal(),
    variable( n ),
    vector_variable( vec_n )
{
    // tell the base class to add the node
    addNode( variable );
    addNode( vector_variable );

}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void NodeTimeSlideUniformProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
NodeTimeSlideUniformProposal* NodeTimeSlideUniformProposal::clone( void ) const
{
    
    return new NodeTimeSlideUniformProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& NodeTimeSlideUniformProposal::getProposalName( void ) const
{
    static std::string name = "NodeTimeSlideUniform";
    
    return name;
}


double NodeTimeSlideUniformProposal::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}


/**
 * Perform the proposal.
 *
 * A Uniform-simplex proposal randomly changes some values of a simplex, although the other values
 * change too because of the renormalization.
 * First, some random indices are drawn. Then, the proposal draws a new somplex
 *   u ~ Uniform(val[index] * alpha)
 * where alpha is the tuning parameter.The new value is set to u.
 * The simplex is then renormalized.
 *
 * \return The hastings ratio.
 */
double NodeTimeSlideUniformProposal::doProposal( void )
{
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;

    // get the current tree either from the single variable or the vector of trees
    Tree *tmp = NULL;
    if ( variable != NULL )
    {
        tmp = &variable->getValue();
    }
    else
    {
        tree_index = floor(rng->uniform01() * vector_variable->getValue().size());
        tmp = &(vector_variable->getValue()[tree_index]);
    }
    Tree& tau = *tmp;

    storedNode = nullptr;

    // pick a random node which is not the root and neithor the direct descendant of the root
    std::vector<TopologyNode*> ok_nodes;
    for(int i=0;i<tau.getNumberOfNodes();i++)
    {
	auto& node = tau.getNode(i);
	if (not node.isRoot() and not node.isTip() and not node.isSampledAncestorTipOrParent())
	    ok_nodes.push_back( & node );
    }

    if (ok_nodes.empty()) return 0.0;

    int index = std::floor(ok_nodes.size() * rng->uniform01());
    TopologyNode* node = ok_nodes[index];
    TopologyNode& parent = node->getParent();

    // we need to work with the times
    double parent_age  = parent.getAge();
    double my_age      = node->getAge();
    double child_Age   = std::max(node->getChild( 0 ).getAge(), node->getChild( 1 ).getAge());

    // now we store all necessary values
    storedNode = node;
    storedAge = my_age;

    // draw new ages and compute the hastings ratio at the same time
    double my_new_age = (parent_age-child_Age) * rng->uniform01() + child_Age;

    // set the age
    tau.getNode(node->getIndex()).setAge( my_new_age );

    return 0.0;
}


/**
 *
 */
void NodeTimeSlideUniformProposal::prepareProposal( void )
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
void NodeTimeSlideUniformProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void NodeTimeSlideUniformProposal::undoProposal( void )
{
    
    // undo the proposal
    if (storedNode)
	storedNode->setAge( storedAge );
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void NodeTimeSlideUniformProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    if ( oldN == variable )
    {
        variable = static_cast<StochasticNode<Tree>* >(newN);
    }
    else if ( oldN == vector_variable )
    {
        vector_variable = static_cast<StochasticNode< RbVector<Tree> >* >( newN );
    }
    
}


void NodeTimeSlideUniformProposal::setProposalTuningParameter(double tp)
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
void NodeTimeSlideUniformProposal::tune( double rate )
{
    
}

