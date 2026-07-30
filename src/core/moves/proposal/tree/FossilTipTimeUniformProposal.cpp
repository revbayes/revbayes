#include <cmath>
#include <iostream>
#include <cstddef>
#include <vector>

#include "DistributionUniform.h"
#include "FossilTipTimeUniformProposal.h"
#include "RbConstants.h"
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
    tip_taxon( t ),
    node_index( 0 ),
    stored_age( 0.0 ),
    failed( false )
{
    // tell the base class to add the node
    addNode( tree );
    addNode( origin );
    addNode( max );
    addNode( min );
    
    if ( tip_taxon == "" )
    {
        use_index = false;
    }
    else
    {
        use_index = true;
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

    bool extended = tree->getDistribution().isExtended();

    failed = false;

    if ( use_index )
    {
        // Always recompute node_index
        node_index = tree->getValue().getTipIndex( tip_taxon );
    }
    else
    {
        std::vector<size_t> tips;
        for (size_t i = 0; i < tau.getNumberOfTips(); ++i)
        {
            TopologyNode* node = &tau.getNode(i);
            // on an extended tree a tip is an extinction, which the taxon declares; its age may
            // already sit at the present, where isFossil() would miss it
            if ( extended == true ? node->getTaxon().isExtinct() : node->isFossil() )
            {
                tips.push_back(i);
            }

        }

        if ( tips.empty() )
        {
            // reject rather than return 0: an unchanged state with ratio 1 is accepted, so a move
            // that can never act reports 100% acceptance and its rate is useless as a diagnostic
            failed = true;
            return RbConstants::Double::neginf;
        }

        // pick a random fossil node
        double u = rng->uniform01();
        node_index = tips[ size_t( std::floor(tips.size() * u) ) ];
    }

    TopologyNode& node = tau.getNode(node_index);

    // an extended tip is bounded by its own occurrences, so it needs no origin
    if ( node.isRoot() == true && origin == NULL && extended == false )
    {
        throw RbException("Attempting to move a root tip, but no origin time provided.");
    }

    // a lone lineage is its own root, so the origin bounds it instead of a parent
    double parent_age   = node.isRoot() ? ( origin != NULL ? origin->getValue() : RbConstants::Double::inf ) : node.getParent().getAge();
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
        TopologyNode& parent = node.getParent();

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

    if ( extended == true )
    {
        // the tip is an extinction, bounded by the present and the youngest occurrence rather than
        // by the parent. A parent-dependent window would need a Hastings term; the density rejects.
        // Occurrence ages are absolute, so this needs the present once it is settable.
        const std::map<TimeInterval, size_t>& occurrences = node.getTaxon().getOccurrences();

        min_age = 0.0;
        max_age = node.getTaxon().getMaxAge();

        // a death lies at or below every occurrence, so the binding bound is the youngest maximum
        for (std::map<TimeInterval, size_t>::const_iterator it = occurrences.begin(); it != occurrences.end(); ++it)
        {
            max_age = fmin(max_age, it->first.getMax());
        }
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
    
    // a proposal that found no movable tip never picked a node, so there is nothing to restore
    if ( failed == true ) return;

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
