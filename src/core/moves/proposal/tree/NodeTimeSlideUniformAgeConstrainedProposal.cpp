#include <cstddef>
#include <cmath>
#include <iostream>
#include <vector>

#include "DistributionUniform.h"
#include "NodeTimeSlideUniformAgeConstrainedProposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbException.h"
#include "Clade.h"
#include "Proposal.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
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
NodeTimeSlideUniformAgeConstrainedProposal::NodeTimeSlideUniformAgeConstrainedProposal( StochasticNode<Tree> *n, const RbVector<Clade> &c ) : Proposal(),
    variable( n ),
    constraints( c )
{
    // tell the base class to add the node
    addNode( variable );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void NodeTimeSlideUniformAgeConstrainedProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
NodeTimeSlideUniformAgeConstrainedProposal* NodeTimeSlideUniformAgeConstrainedProposal::clone( void ) const
{
    
    return new NodeTimeSlideUniformAgeConstrainedProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& NodeTimeSlideUniformAgeConstrainedProposal::getProposalName( void ) const
{
    static std::string name = "NodeTimeSlideUniformAgeConstrained";
    
    return name;
}


double NodeTimeSlideUniformAgeConstrainedProposal::getProposalTuningParameter( void ) const
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
double NodeTimeSlideUniformAgeConstrainedProposal::doProposal( void )
{
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    Tree& tau = variable->getValue();
    
    if (tau.getNumberOfTips() <= 2)
    {
        return 0.0;
    }
    
    // pick a random node which is not the root and neithor the direct descendant of the root
    std::vector<TopologyNode*> nodes;
    double min_age = -1;
    double max_age = -1;
    double our_age = -1;
    for (size_t i=0; i<constraints.size(); ++i)
    {
        TopologyNode *node = &tau.getMrca( constraints[i] );
        
        TopologyNode& parent = node->getParent();
        
        // we need to work with the times
        double parent_age  = parent.getAge();
        double my_age      = node->getAge();
        double child_Age   = node->getChild( 0 ).getAge();
        if ( child_Age < node->getChild( 1 ).getAge())
        {
            child_Age = node->getChild( 1 ).getAge();
        }
        
        if ( min_age == -1 || min_age < child_Age )
        {
            min_age = child_Age;
        }
        
        if ( max_age == -1 || max_age > parent_age )
        {
            max_age = parent_age;
        }
        
        if ( our_age == -1 )
        {
            our_age = my_age;
        }
        else if ( fabs( our_age - my_age ) > 1E-10 )
        {
            throw RbException("Problem in joint node age update because node ages are not equal.");
        }
        
        
        nodes.push_back( node );
    }
    
    
    // now we store all necessary values
    stored_age = our_age;
    
    // draw new ages and compute the hastings ratio at the same time
    double our_new_age = (max_age-min_age) * rng->uniform01() + min_age;
    
    
    // set the age
    for (size_t i=0; i<constraints.size(); ++i)
    {
        TopologyNode *node = nodes[i];
        node->setAge( our_new_age );
    }
    
    
    return 0.0;
    
}


/**
 *
 */
void NodeTimeSlideUniformAgeConstrainedProposal::prepareProposal( void )
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
void NodeTimeSlideUniformAgeConstrainedProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void NodeTimeSlideUniformAgeConstrainedProposal::undoProposal( void )
{
    
    Tree& tau = variable->getValue();
    
    // undo the proposal
    for (size_t i=0; i<constraints.size(); ++i)
    {
        TopologyNode *node = &tau.getMrca( constraints[i] );
        node->setAge( stored_age );
    }
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void NodeTimeSlideUniformAgeConstrainedProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    variable = static_cast<StochasticNode<Tree>* >(newN) ;
    
}


void NodeTimeSlideUniformAgeConstrainedProposal::setProposalTuningParameter(double tp)
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
void NodeTimeSlideUniformAgeConstrainedProposal::tune( double rate )
{
    
}

