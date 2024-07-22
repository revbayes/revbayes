#include <cmath>
#include <iostream>
#include <cstddef>
#include <vector>

#include "DistributionUniform.h"
#include "FossilTipTimeSlideUniformProposal.h"
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
FossilTipTimeSlideUniformProposal::FossilTipTimeSlideUniformProposal( StochasticNode<Tree> *n, TypedDagNode<double> *o, TypedDagNode<double> *ma, TypedDagNode<double> *mi, const std::string& t, double l, double r ) : Proposal(r),
    tree( n ),
    origin( o ),
    max( ma ),
    min( mi ),
    tip_taxon( t ),
    lambda( l )
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
void FossilTipTimeSlideUniformProposal::cleanProposal( void )
{
    // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
FossilTipTimeSlideUniformProposal* FossilTipTimeSlideUniformProposal::clone( void ) const
{
    
    return new FossilTipTimeSlideUniformProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& FossilTipTimeSlideUniformProposal::getProposalName( void ) const
{
    static std::string name = "FossilTipTimeSlideUniform";
    
    return name;
}


double FossilTipTimeSlideUniformProposal::getProposalTuningParameter( void ) const
{
    // get the scaling parameter for the proposal
    return lambda;
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
double FossilTipTimeSlideUniformProposal::doProposal( void )
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
    
    TopologyNode* node = &tau.getNode(node_index);

    TopologyNode& parent = node->getParent();

    // we need to work with the times
    double parent_age   = parent.getAge();
    double my_age       = node->getAge();    
    double min_age      = 0;
    double max_age      = parent_age;
    
    // adjust min and max age, either given taxon data or given provided ages
    if ( min == NULL )
    {
        // adjust min age given taxon data
        Taxon& taxon = node->getTaxon();
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
        Taxon& taxon = node->getTaxon();
        double taxon_max_age = taxon.getMaxAge();
        if ( taxon_max_age < max_age )
        {
            max_age = taxon_max_age;
        }
    }
    else
    {
        // adjust max age given provided variable
        double provided_max_age = max->getValue();
        if ( provided_max_age < max_age )
        {
            max_age = provided_max_age;
        }
    }

    if ( node->isSampledAncestorTip() == true )
    {
        TopologyNode *sibling = &parent.getChild( 0 );
        if ( sibling == node )
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
    }
    
    // now we store all necessary values
    stored_age = my_age;
    
    double size = max_age - min_age;
    
    double u      = rng->uniform01();
    double delta  = ( lambda * ( u - 0.5 ) );
    
    if ( fabs(delta) > 2.0*size )
    {
        delta -= floor(delta / (2.0*size)) * (2.0*size);
    }
    double new_age = my_age + delta;
    
    /* reflect the new value */
    do {
        if ( new_age < min_age )
        {
            new_age = 2.0 * min_age - new_age;
        }
        else if ( new_age > max_age )
        {
            new_age = 2.0 * max_age - new_age;
        }
    } while ( new_age < min_age || new_age > max_age );
    
    // set the age
    node->setAge( new_age );
    
    
    // this is a symmetric proposal so the hasting ratio is 0.0
    return 0.0;
}


/**
 *
 */
void FossilTipTimeSlideUniformProposal::prepareProposal( void )
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
void FossilTipTimeSlideUniformProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
    o << "delta = ";
    if (name_only == false)
    {
        o << lambda;
    }
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void FossilTipTimeSlideUniformProposal::undoProposal( void )
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
void FossilTipTimeSlideUniformProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    if (oldN == tree)
    {
        tree = static_cast<StochasticNode<Tree>* >(newN) ;
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


void FossilTipTimeSlideUniformProposal::setProposalTuningParameter(double tp)
{
    lambda = tp;
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
void FossilTipTimeSlideUniformProposal::tune( double rate )
{
    
    double p = this->targetAcceptanceRate;
    if ( rate > p )
    {
        lambda *= (1.0 + ((rate-p)/(1.0 - p)));
    }
    else
    {
        lambda /= (2.0 - rate/p);
    }
    
    lambda = fmin(10000, lambda);
    
}
