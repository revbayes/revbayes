#include <cmath>
#include <iostream>
#include <cstddef>
#include <vector>

#include "DistributionUniform.h"
#include "FossilTipTimeUniformProposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbException.h"
#include "TypedDagNode.h"
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
FossilTipTimeUniformProposal::FossilTipTimeUniformProposal( StochasticNode<Tree> *n, TypedDagNode<double> *o, TypedDagNode<double> *ma, TypedDagNode<double> *mi, const std::string& t ) : Proposal(),
    tree( n ),
    origin( o ),
    max( ma ),
    min( mi ),
    tip_taxon( t )
{
    // tell the base class to add the node
    addNode( tree );
    addNode( origin );
    addNode( max );
    addNode( min );
    
    if ( tip_taxon == "" )
    {
        use_index = false;
        node_index = -1;
    }
    else
    {
        use_index = true;
        node_index = tree->getValue().getTipIndex( tip_taxon );
    }
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void FossilTipTimeUniformProposal::cleanProposal( void )
{
    // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
FossilTipTimeUniformProposal* FossilTipTimeUniformProposal::clone( void ) const
{
    
    return new FossilTipTimeUniformProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& FossilTipTimeUniformProposal::getProposalName( void ) const
{
    static std::string name = "FossilTipTimeUniform";
    
    return name;
}


double FossilTipTimeUniformProposal::getProposalTuningParameter( void ) const
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
double FossilTipTimeUniformProposal::doProposal( void )
{
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    Tree& tau = tree->getValue();
    
    if ( use_index == false )
    {
        std::vector<size_t> tips;
        for (size_t i = 0; i < tau.getNumberOfTips(); ++i)
        {
            TopologyNode* node = &tau.getNode(i);
            if ( node->isFossil() )
            {
                tips.push_back(i);
            }

        }

        if ( tips.empty() )
        {
            failed = true;
            return 0;
        }

        // pick a random fossil node
        double u = rng->uniform01();
        node_index = tips[ size_t( std::floor(tips.size() * u) ) ];
    }

    TopologyNode& node = tau.getNode(node_index);
    TopologyNode& parent = node.getParent();

    // we need to work with the times
    double parent_age   = parent.getAge();
    double my_age       = node.getAge();
    double min_age      = 0;
    double max_age;
    
    // adjust min and max age, either given taxon data or given provided ages
    if ( min == NULL )
    {
        // adjust min age given taxon data
        Taxon& taxon = node.getTaxon();
        min_age = taxon.getMinAge();
    }
    else
    {
        // adjust min age given provided age
        min_age = min->getValue();
    }
    if ( max == NULL )
    {
        // adjust max age given taxon data
        Taxon& taxon = node.getTaxon();
        double taxon_max_age = taxon.getMaxAge();
        max_age = taxon_max_age;
    }
    else
    {
        // adjust max age given provided variable
        double provided_max_age = max->getValue();
        max_age = provided_max_age;
    }

    if ( node.isSampledAncestorTip() == true )
    {
        TopologyNode *sibling = &parent.getChild( 0 );
        if ( sibling == &node )
        {
            sibling = &parent.getChild( 1 );
        }

        double sib_age = sibling->getAge();
        min_age = fmax(min_age, sib_age);

        if ( parent.isRoot() )
        {
            if (origin == NULL)
            {
                throw RbException("Attempting to move root sampled ancestor, but no origin time provided.");
            }
            
            double origin_age = origin->getValue();
            
            // set the max age either to the boundary or the parent max age
            max_age = fmin(max_age, origin_age);
        }
        else
        {
            TopologyNode& grandParent = parent.getParent();

            double grandparent_age = grandParent.getAge();
            
            // set the max age either to the boundary or the parent max age
            max_age = fmin(max_age, grandparent_age);
        }
    } else {
        max_age = fmin(max_age, parent_age);
    }
    
    assert(max_age >= min_age); //sanity check

    // now we store all necessary values
    stored_age = my_age;
    
    // draw new ages and compute the hastings ratio at the same time
    double my_new_age = min_age + (max_age - min_age) * rng->uniform01();
    
    // set the age
    node.setAge( my_new_age );

    return 0.0;
}


/**
 *
 */
void FossilTipTimeUniformProposal::prepareProposal( void )
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
void FossilTipTimeUniformProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void FossilTipTimeUniformProposal::undoProposal( void )
{
    
    // undo the proposal
    Tree& tau = tree->getValue();
    TopologyNode* node = &tau.getNode(node_index);
    node->setAge( stored_age );
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void FossilTipTimeUniformProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    if (oldN == tree)
    {
        tree = static_cast<StochasticNode<Tree>* >(newN);
        if ( tip_taxon != "" )
        {
            node_index = tree->getValue().getTipIndex( tip_taxon );
        }
    }
    else if (oldN == origin)
    {
        origin = static_cast<TypedDagNode<double>* >(newN) ;
    }
    else if (oldN == max)
    {
        max = static_cast<TypedDagNode<double>* >(newN) ;
    }
    else if (oldN == min)
    {
        min = static_cast<TypedDagNode<double>* >(newN) ;
    }
    
}


void FossilTipTimeUniformProposal::setProposalTuningParameter(double tp)
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
void FossilTipTimeUniformProposal::tune( double rate )
{
    
}
