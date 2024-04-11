#include <cmath>
#include <iostream>
#include <cstdlib>
#include <vector>

#include "DistributionUniform.h"
#include "AddRemoveTipProposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
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
AddRemoveTipProposal::AddRemoveTipProposal( StochasticNode<Tree> *n, bool exa, bool exi, bool sa) : Proposal(),
    tau( n ),
    storedSibling(NULL),
    storedTip(NULL),
    removed(false),
    extant(exa),
    extinct(exi),
    failed(false),
    sampled_ancestors(sa)
{
    // tell the base class to add the node
    addNode( tau );

}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void AddRemoveTipProposal::cleanProposal( void )
{
    if (removed)
    {
        removed = false;
        delete &storedTip->getParent();
    }

    failed = false;
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
AddRemoveTipProposal* AddRemoveTipProposal::clone( void ) const
{
    
    return new AddRemoveTipProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& AddRemoveTipProposal::getProposalName( void ) const
{
    static std::string name = "AddRemoveTip";
    
    return name;
}


double AddRemoveTipProposal::getProposalTuningParameter( void ) const
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
double AddRemoveTipProposal::doProposal( void )
{
    failed = false;

    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;

    Tree &t = tau->getValue();

    std::vector<TopologyNode*> tips;
    std::vector<TopologyNode*> siblings;
    for (size_t i = 0; i < t.getNumberOfNodes(); ++i)
    {
        TopologyNode* node = &t.getNode(i);

        if ( node->isTip() == true && node->getParent().isRoot() == false )
        {
            if ( ( extinct == true && node->isFossil() == true && node->isSampledAncestorTip() == sampled_ancestors ) ||
                (  extant == true && node->isFossil() == false ) )
            {
                if( ( node->isSampledAncestorTip() == false && node->getParent().isSampledAncestorTipOrParent() == true ) == false )
                    tips.push_back(node);
            }
        }

        if ( node->isRoot() == false && node->isSampledAncestorTip() == false )
        {
            siblings.push_back(node);
        }

    }

    TopologyNode* node;

    double u = rng->uniform01();

    double hr = 0;
    double jacobian = 0;

    bool both = extinct and extant;

    // pick a random tip node to remove
    if (u < 0.5)
    {
        if (tips.empty())
        {
            failed = true;
            return 0.0;
        }

        u = rng->uniform01();
        node = tips[ size_t(u*tips.size()) ];

        hr = log( tips.size() / double(siblings.size()-2) ) - both*RbConstants::LN2;
        jacobian = removeTip( node );
    }
    // or just pick a random node above which to add a tip
    else
    {
        u = rng->uniform01();
        node = siblings[ size_t(u*siblings.size()) ];

        hr = log( siblings.size() / double(tips.size()+1) ) + both*RbConstants::LN2;
        jacobian = addTip( node );
    }

    return hr + jacobian;
    
}



/*
 * Add a fossil, either sampled ancestor or not
 */
double AddRemoveTipProposal::addTip(TopologyNode *n)
{
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;

    Tree &t = tau->getValue();

    storedSibling = n;

    std::stringstream name;
    name << "t" << std::rand();
    storedTip = new TopologyNode(name.str());
    storedTip->setAge(0.0);

    TopologyNode* new_parent_node = new TopologyNode();
    new_parent_node->addChild(storedTip);
    storedTip->setParent(new_parent_node);

    TopologyNode* grandParent = &storedSibling->getParent();
    grandParent->removeChild(storedSibling);
    grandParent->addChild(new_parent_node);
    new_parent_node->setParent(grandParent);
    new_parent_node->addChild(storedSibling);
    storedSibling->setParent(new_parent_node);

    t.setRoot(&t.getRoot(), true);
    
    // determine upper/lower bound
    double min_age = storedSibling->getAge();
    double max_age = grandParent->getAge();

    // set the new age for the parent node
    double u = rng->uniform01();
    new_parent_node->setAge( (max_age - min_age) * u + min_age );
    
    double lnJacobian = log(max_age - min_age);

    // not yet sure why this is needed, but it is
    // unless lnProbTreeShape is incorrect for all birth death processes
    double hr = 0;

    if ( extinct and extant )
    {
        u = rng->uniform01();
    }
    else
    {
        u = extant;
    }

    if ( u < 0.5 )
    {
        if ( sampled_ancestors == false )
        {
            double v = rng->uniform01();

            lnJacobian += log(new_parent_node->getAge());

            // draw the new age for the fossil node
            double new_tip_age = new_parent_node->getAge() * v;

            storedTip->setAge(new_tip_age);
        }
        else
        {
            storedTip->setAge( new_parent_node->getAge() );
            storedTip->setSampledAncestor(true);
        }

        //hr = log( t.getNumberOfExtinctTips() );
    }
    else
    {
        hr = log( t.getNumberOfExtantTips() );
    }

    return lnJacobian + hr;
}



/*
 * Remove a fossil
 */
double AddRemoveTipProposal::removeTip(TopologyNode *n)
{
    removed = true;

    Tree &t = tau->getValue();

    storedTip = n;

    TopologyNode* parent = &storedTip->getParent();
    TopologyNode* grandParent = &parent->getParent();

    storedSibling = &parent->getChild( 0 );
    if ( storedSibling == n )
    {
        storedSibling = &parent->getChild( 1 );
    }

    parent->removeChild(storedSibling);
    grandParent->removeChild(parent);
    grandParent->addChild(storedSibling);
    storedSibling->setParent(grandParent);

    parent->setParent(NULL);

    t.setRoot(&t.getRoot(), true);


    double max_age = grandParent->getAge();
    double min_age = storedSibling->getAge();

    double lnJacobian = - log(max_age - min_age);

    // not yet sure why this is needed, but it is
    // unless lnProbTreeShape is incorrect for all birth death processes
    double hr = 0;

    if ( storedTip->getAge() > 0.0 )
    {
        if ( sampled_ancestors == false )
        {
            lnJacobian -= log(parent->getAge());
        }

        //hr = - log( t.getNumberOfExtinctTips() + 1 );
    }
    else
    {
        hr = - log( t.getNumberOfExtantTips() + 1 );
    }

    return lnJacobian + hr;
}



/**
 *
 */
void AddRemoveTipProposal::prepareProposal( void )
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
void AddRemoveTipProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void AddRemoveTipProposal::undoProposal( void )
{
    if ( failed == true )
        return;

    Tree &t = tau->getValue();

    TopologyNode *parent = &storedTip->getParent();

    // undo the add proposal
    if ( removed == false )
    {
        TopologyNode* grandParent = &parent->getParent();

        parent->removeChild(storedSibling);
        grandParent->removeChild(parent);
        grandParent->addChild(storedSibling);
        storedSibling->setParent(grandParent);

        parent->setParent(NULL);

        t.setRoot(&t.getRoot(), true);

        delete parent;
    }
    // undo the remove proposal
    else
    {
        removed = false;

        TopologyNode* grandParent = &storedSibling->getParent();

        grandParent->removeChild(storedSibling);
        grandParent->addChild(parent);
        parent->setParent(grandParent);
        parent->addChild(storedSibling);
        storedSibling->setParent(parent);

        t.setRoot(&t.getRoot(), true);
    }

}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void AddRemoveTipProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    if ( oldN == tau )
    {
        tau = static_cast<StochasticNode<Tree>* >(newN) ;
    }
    
}


void AddRemoveTipProposal::setProposalTuningParameter(double tp)
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
void AddRemoveTipProposal::tune( double rate )
{
    
}

