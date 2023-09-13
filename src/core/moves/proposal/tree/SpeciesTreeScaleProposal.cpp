#include "SpeciesTreeScaleProposal.h"

#include <cstddef>
#include <cmath>
#include <iostream>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "TreeUtilities.h"
#include "Cloneable.h"
#include "DagNode.h"
#include "StochasticNode.h"
#include "TopologyNode.h"
#include "Tree.h"

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
SpeciesTreeScaleProposal::SpeciesTreeScaleProposal( StochasticNode<Tree> *sp, StochasticNode<double> *r, double d ) : Proposal(),
    speciesTree( sp ),
    geneTrees( ),
    rootAge( r ),
    delta( d )
{

    // tell the base class to add the node
    addNode( speciesTree );
    addNode( rootAge );

}


SpeciesTreeScaleProposal::~SpeciesTreeScaleProposal( void )
{

}


/**
 * Add a new DAG node holding a gene tree on which this move operates on.
 *
 */
void SpeciesTreeScaleProposal::addGeneTree(StochasticNode<Tree> *gt)
{
    // check if this node isn't already in our list
    bool exists = false;
    for (size_t i=0; i < geneTrees.size(); ++i)
    {
        if ( geneTrees[i] == gt )
        {
            exists = true;
            break;
        }
    }

    // only add this variable if it doesn't exist in our list already
    if ( exists == false )
    {
        geneTrees.push_back( gt );
        addNode( gt );
    }

}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void SpeciesTreeScaleProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
SpeciesTreeScaleProposal* SpeciesTreeScaleProposal::clone( void ) const
{

    return new SpeciesTreeScaleProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& SpeciesTreeScaleProposal::getProposalName( void ) const
{
    static std::string name = "SpeciesTreeScale";

    return name;
}


double SpeciesTreeScaleProposal::getProposalTuningParameter( void ) const
{
    return delta;
}


/**
 * Perform the proposal.
 *
 * \return The hastings ratio.
 */
double SpeciesTreeScaleProposal::doProposal( void )
{

    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;

    Tree& tau = speciesTree->getValue();

    // get the root of the species tree
    TopologyNode& root = tau.getRoot();

    // we need to work with the times
    double my_age      = root.getAge();

    // now we store all necessary values
    storedAge = my_age;

    // draw new ages and compute the hastings ratio at the same time
    double u = rng->uniform01();
    double scaling_factor = std::exp( delta * ( u - 0.5 ) );

    // Sebastian: This is for debugging to test if the proposal's acceptance rate is 1.0 as it should be!
//    scaling_factor = 1.0;

    double my_new_age = my_age * scaling_factor;

    if ( rootAge != NULL )
    {
        rootAge->setValue( new double(my_new_age) );
    }

    size_t num_nodes = tau.getNumberOfInteriorNodes();

    for ( size_t i=0; i<geneTrees.size(); ++i )
    {
        // get the i-th gene tree
        Tree& gene_tree = geneTrees[i]->getValue();

        // rescale the subtree of this gene tree
        TreeUtilities::rescaleSubtree(gene_tree.getRoot(), scaling_factor );

        num_nodes += gene_tree.getNumberOfInteriorNodes();

        // Sebastian: This is only for debugging. It makes the code slower. Hopefully it is not necessary anymore.
//        geneTrees[i]->touch( true );

    }

    // Sebastian: We need to work on a mechanism to make these proposal safe for non-ultrametric trees!
    //    if (min_age != 0.0)
    //    {
    //        for (size_t i = 0; i < tau.getNumberOfTips(); i++)
    //        {
    //            if (tau.getNode(i).getAge() < 0.0)
    //            {
    //                return RbConstants::Double::neginf;
    //            }
    //        }
    //    }


    // rescale the subtree of the species tree
    TreeUtilities::rescaleSubtree(root, scaling_factor );

    // compute the Hastings ratio
    double lnHastingsratio = log( scaling_factor ) * num_nodes;


    return lnHastingsratio;
}


/**
 *
 */
void SpeciesTreeScaleProposal::prepareProposal( void )
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
void SpeciesTreeScaleProposal::printParameterSummary(std::ostream &o, bool name_only) const
{

    o << "delta = ";
    if (name_only == false)
    {
        o << delta;
    }

}


/**
 * Remove a DAG node holding a gene tree on which this move operates on.
 *
 */
void SpeciesTreeScaleProposal::removeGeneTree(StochasticNode<Tree> *gt)
{
    // remove it from our list
    for (size_t i=0; i < geneTrees.size(); ++i)
    {
        if ( geneTrees[i] == gt )
        {
            geneTrees.erase( geneTrees.begin() + i );
            --i;
        }
    }

}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void SpeciesTreeScaleProposal::undoProposal( void )
{
    // undo the proposal

    Tree& tau = speciesTree->getValue();

    TopologyNode& node = tau.getRoot();

    double sf = storedAge / node.getAge();

    // re-scale
    TreeUtilities::rescaleSubtree(node, sf );

    if ( rootAge != NULL )
    {
        rootAge->setValue( new double(storedAge) );
    }

    for ( size_t i=0; i<geneTrees.size(); ++i )
    {
        // get the i-th gene tree
        Tree& gene_tree = geneTrees[i]->getValue();

        // rescale the subtree of this gene tree
        //TreeUtilities::rescaleTree(&gene_tree, &gene_tree.getRoot(), sf );
        TreeUtilities::rescaleSubtree(gene_tree.getRoot(), sf );

    }

}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void SpeciesTreeScaleProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{

    if ( oldN == rootAge )
    {
        rootAge = static_cast<StochasticNode<double>* >(newN) ;
    }
    else if ( oldN == speciesTree )
    {
        speciesTree = static_cast<StochasticNode<Tree>* >(newN) ;
    }
    else
    {
        size_t num_found_trees = 0;
        for ( size_t i=0; i<geneTrees.size(); ++i )
        {
            if ( oldN == geneTrees[i] )
            {
                geneTrees[i] = static_cast<StochasticNode<Tree>* >(newN) ;
                ++num_found_trees;
            }
        }

        if ( num_found_trees != 1 )
        {
            std::cerr << "Found " << num_found_trees << " DAG nodes matching to node \"" << oldN->getName() << "\".";
        }

    }

}


void SpeciesTreeScaleProposal::setProposalTuningParameter(double tp)
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
void SpeciesTreeScaleProposal::tune( double rate )
{

    if ( rate > 0.234 )
    {
        delta *= (1.0 + ((rate-0.234)/0.766) );
    }
    else
    {
        delta /= (2.0 - rate/0.234 );
    }

}
