#ifndef FixedNodeheightPruneAndRegraftCharacterHistoryProposal_H
#define FixedNodeheightPruneAndRegraftCharacterHistoryProposal_H

#include "BranchHistoryDiscrete.h"
#include "HomologousDiscreteCharacterData.h"
#include "PathRejectionSampleProposal.h"
#include "Proposal.h"
#include "NodeRejectionSampleProposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RateGenerator.h"
#include "RateGeneratorSequence.h"
#include "RbConstants.h"
#include "RbException.h"
#include "StochasticNode.h"
#include "Tree.h"
#include "TypedDagNode.h"

#include <cmath>
#include <iostream>
#include <set>
#include <string>

namespace RevBayesCore {
    
    /**
     * The fixed node-height prune-and-regraft operator.
     *
     * A fixed node-height prune-and-regraft proposal is a SPR (subtree prune-and-regraft) proposal on rooted trees without changing the node age.
     * That is, we pick a random node which is not the root.
     * Then, we prune this node and try to attach it anywhere else in the tree at this node age.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2012-07-12, version 1.0
     *
     */
    
    template<class charType>
    class FixedNodeheightPruneAndRegraftCharacterHistoryProposal : public Proposal {
        
    public:
        FixedNodeheightPruneAndRegraftCharacterHistoryProposal( StochasticNode<Tree> *t, StochasticNode<AbstractHomologousDiscreteCharacterData>* c );                                   //!<  constructor
        
        // Basic utility functions
        void                                                        cleanProposal(void);                                        //!< Clean up proposal
        FixedNodeheightPruneAndRegraftCharacterHistoryProposal*     clone(void) const;                                          //!< Clone object
        double                                                      doProposal(void);                                           //!< Perform proposal
        const std::string&                                          getProposalName(void) const;                                //!< Get the name of the proposal for summary printing
        double                                                      getProposalTuningParameter(void) const;
        void                                                        prepareProposal(void);                                      //!< Prepare the proposal
        void                                                        printParameterSummary(std::ostream &o, bool name_only) const;               //!< Print the parameter summary
        void                                                        sampleNodeCharacters(TopologyNode* node);                       //!< Sample the characters at the node
        void                                                        setRateGenerator(const TypedDagNode<RateGenerator> *d);
        void                                                        setRateGenerator(const TypedDagNode<RateGeneratorSequence> *d);
        void                                                        setProposalTuningParameter(double tp);
        void                                                        tune(double r);                                             //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                                        undoProposal(void);                                         //!< Reject the proposal
        
    protected:
        
        void                                                        swapNodeInternal(DagNode *oldN, DagNode *newN);             //!< Swap the DAG nodes on which the Proposal is working on
        
    private:
        // private helper methods
        void                                                        findNewBrothers(std::vector<TopologyNode*> &b, TopologyNode &p, TopologyNode *n);

        // parameters
        StochasticNode<Tree>*                                       tree;                                                           //!< The tree the Proposal is working on
        StochasticNode<AbstractHomologousDiscreteCharacterData>*    ctmc;                                                           //!< The character history the Proposal is working on
        const TypedDagNode<RateGenerator>*                          rate_generator;
        const TypedDagNode<RateGeneratorSequence>*                  rate_generator_sequence;

        // dimensions
        size_t                                                      num_states;

        // proposals
        PathRejectionSampleProposal<charType>*                      node_proposal_one;
        PathRejectionSampleProposal<charType>*                      left_proposal_one;
        PathRejectionSampleProposal<charType>*                      right_proposal_one;
        PathRejectionSampleProposal<charType>*                      node_proposal_two;
        PathRejectionSampleProposal<charType>*                      left_proposal_two;
        PathRejectionSampleProposal<charType>*                      right_proposal_two;
        PathRejectionSampleProposal<charType>*                      left_proposal_three;
        PathRejectionSampleProposal<charType>*                      right_proposal_three;

        // transition probability matrices
        TransitionProbabilityMatrix                                 node_tp_matrix;
        TransitionProbabilityMatrix                                 left_tp_matrix;
        TransitionProbabilityMatrix                                 right_tp_matrix;

        // stored objects to undo proposal
        bool                                                        failed;
        TopologyNode*                                               stored_brother;
        TopologyNode*                                               stored_new_brother;

        std::vector<size_t>                                         stored_node_states_one;
        std::vector<size_t>                                         stored_node_states_two;
        std::vector<size_t>                                         stored_node_states_three;
        
    };
    
}

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
template<class charType>
RevBayesCore::FixedNodeheightPruneAndRegraftCharacterHistoryProposal<charType>::FixedNodeheightPruneAndRegraftCharacterHistoryProposal( StochasticNode<Tree> *t, StochasticNode<AbstractHomologousDiscreteCharacterData>* c ) : Proposal(),
    tree( t ),
    ctmc( c ),
    rate_generator( NULL ),
    rate_generator_sequence( NULL ),
    num_states( 2 ),
    node_tp_matrix( 2 ),
    left_tp_matrix( 2 ),
    right_tp_matrix( 2 )
{

    // tell the base class to add the nodes
    addNode( tree );
    addNode( ctmc );
    
    // make the proposals
    node_proposal_one  = new PathRejectionSampleProposal<charType>( c );
    left_proposal_one  = new PathRejectionSampleProposal<charType>( c );
    right_proposal_one = new PathRejectionSampleProposal<charType>( c );
    
    node_proposal_two  = new PathRejectionSampleProposal<charType>( c );
    left_proposal_two  = new PathRejectionSampleProposal<charType>( c );
    right_proposal_two = new PathRejectionSampleProposal<charType>( c );
    
    left_proposal_three  = new PathRejectionSampleProposal<charType>( c );
    right_proposal_three = new PathRejectionSampleProposal<charType>( c );

}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
template<class charType>
void RevBayesCore::FixedNodeheightPruneAndRegraftCharacterHistoryProposal<charType>::cleanProposal( void )
{
    
    node_proposal_one->cleanProposal();
    left_proposal_one->cleanProposal();
    right_proposal_one->cleanProposal();
    
    node_proposal_two->cleanProposal();
    left_proposal_two->cleanProposal();
    right_proposal_two->cleanProposal();

    left_proposal_three->cleanProposal();
    right_proposal_three->cleanProposal();

    stored_node_states_one.clear();
    stored_node_states_two.clear();
    stored_node_states_three.clear();
    
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
template<class charType>
RevBayesCore::FixedNodeheightPruneAndRegraftCharacterHistoryProposal<charType>* RevBayesCore::FixedNodeheightPruneAndRegraftCharacterHistoryProposal<charType>::clone( void ) const
{
    
    return new FixedNodeheightPruneAndRegraftCharacterHistoryProposal( *this );
}


template<class charType>
void RevBayesCore::FixedNodeheightPruneAndRegraftCharacterHistoryProposal<charType>::findNewBrothers(std::vector<TopologyNode *> &b, TopologyNode &p, TopologyNode *n) {
    // security check that I'm not a tip
    if (!n->isTip() && &p != n)
    {
        // check the first child
        TopologyNode* child = &n->getChild( 0 );
        if ( child->getAge() < p.getAge() )
        {
            if ( child->isSampledAncestorTip() == false ) b.push_back( child );
        }
        else
        {
            findNewBrothers(b, p, child);
        }
        
        // check the second child
        child = &n->getChild( 1 );
        if ( child->getAge() < p.getAge() )
        {
            if ( child->isSampledAncestorTip() == false ) b.push_back( child );
        }
        else
        {
            findNewBrothers(b, p, child);
        }
    }
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
template<class charType>
const std::string& RevBayesCore::FixedNodeheightPruneAndRegraftCharacterHistoryProposal<charType>::getProposalName( void ) const
{
    static std::string name = "FNPRCharacterHistory";
    
    return name;
}


template<class charType>
double RevBayesCore::FixedNodeheightPruneAndRegraftCharacterHistoryProposal<charType>::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}


/**
 * Perform the proposal.
 *
 * A Beta-simplex proposal randomly changes some values of a simplex, although the other values
 * change too because of the renormalization.
 * First, some random indices are drawn. Then, the proposal draws a new somplex
 *   u ~ Beta(val[index] * alpha)
 * where alpha is the tuning parameter.The new value is set to u.
 * The simplex is then renormalized.
 *
 * \return The hastings ratio.
 */
template<class charType>
double RevBayesCore::FixedNodeheightPruneAndRegraftCharacterHistoryProposal<charType>::doProposal( void )
{
    
    // reset flags
    failed = false;
    
    // Get random number generator
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    Tree& tau = tree->getValue();
    
    // pick a random node which is not the root and neithor the direct descendant of the root
    TopologyNode* node;
    do {
        double u = rng->uniform01();
        size_t index = size_t( std::floor(tau.getNumberOfNodes() * u) );
        node = &tau.getNode(index);
    } while ( node->isRoot() || node->getParent().isRoot() || node->getParent().getChild(0).isSampledAncestorTip() || node->getParent().getChild(1).isSampledAncestorTip() );
    
    TopologyNode* parent        = &node->getParent();
    TopologyNode* grandparent   = &parent->getParent();
    TopologyNode* brother       = &parent->getChild( 0 );
    // check if we got the correct child
    if ( brother == node )
    {
        brother = &parent->getChild( 1 );
    }
    
    // collect the possible reattachement points
    std::vector<TopologyNode*> new_brothers;
    findNewBrothers(new_brothers, *parent, &tau.getRoot());
    
    // we only need to propose a new tree if there are any other re-attachement points
    if ( new_brothers.size() < 1 )
    {
        failed = true;
        return RbConstants::Double::neginf;
    }

    size_t index = size_t(rng->uniform01() * new_brothers.size());
    TopologyNode* new_brother = new_brothers[index];
    TopologyNode* new_grandparent = &new_brother->getParent();
    
    // step 0: retrieve the character history
    TreeHistoryCtmc<charType>* p = dynamic_cast< TreeHistoryCtmc<charType>* >( &ctmc->getDistribution() );
    if ( p == NULL )
    {
        throw RbException("Failed cast.");
    }
    size_t num_sites = p->getNumberOfSites();

//    std::cout << std::endl;
//    std::cout << "Initial probability: " << p->computeLnProbability() << std::endl;
    
    // step 1: store all the necessary values
    stored_brother     = brother;
    stored_new_brother = new_brother;

    stored_node_states_one.clear();
    const std::vector<CharacterEvent*>& node_states_one = p->getHistory( *grandparent ).getChildCharacters();
    stored_node_states_one.resize(num_sites, 0);
//    std::cout << "Initial states one: ";
    for (size_t site_index = 0; site_index < num_sites; ++site_index)
    {
        size_t s = static_cast<CharacterEventDiscrete*>(node_states_one[site_index])->getState();
        stored_node_states_one[site_index] = s;
//        std::cout << s << " ";
    }
//    std::cout << std::endl;
    
    stored_node_states_two.clear();
    const std::vector<CharacterEvent*>& node_states_two = p->getHistory( *new_grandparent ).getChildCharacters();
    stored_node_states_two.resize(num_sites, 0);
//    std::cout << "Initial states two: ";
    for (size_t site_index = 0; site_index < num_sites; ++site_index)
    {
        size_t s = static_cast<CharacterEventDiscrete*>(node_states_two[site_index])->getState();
        stored_node_states_two[site_index] = s;
//        std::cout << s << " ";
    }
//    std::cout << std::endl;

    stored_node_states_three.clear();
    const std::vector<CharacterEvent*>& node_states_three = p->getHistory( *parent ).getChildCharacters();
    stored_node_states_three.resize(num_sites, 0);
//    std::cout << "Initial states three: ";
    for (size_t site_index = 0; site_index < num_sites; ++site_index)
    {
        size_t s = static_cast<CharacterEventDiscrete*>(node_states_three[site_index])->getState();
        stored_node_states_three[site_index] = s;
//        std::cout << s << " ";
    }
//    std::cout << std::endl;

    // step 2: exchange the two nodes
    grandparent->removeChild( parent );
    parent->removeChild( brother );
    grandparent->addChild( brother );               // need to sample a new state for grandparent
    brother->setParent( grandparent );
    
    new_grandparent->removeChild( new_brother );
    new_grandparent->addChild( parent );            // need to sample a new state for new_grandparent
    parent->setParent( new_grandparent );
    parent->addChild( new_brother );                // need to sample a new state for parent
    new_brother->setParent( parent );
    
    // step 3: set up the path proposals, making sure
    // not to set two proposals for the same branch
    node_proposal_one->assignNode( grandparent );
    node_proposal_one->prepareProposal();
    left_proposal_one->assignNode( &grandparent->getChild(0) );
    left_proposal_one->prepareProposal();
    right_proposal_one->assignNode( &grandparent->getChild(1) );
    right_proposal_one->prepareProposal();

    if ( &new_grandparent != &grandparent && new_grandparent != &grandparent->getChild(0) && new_grandparent != &grandparent->getChild(1) )
    {
        node_proposal_two->assignNode( new_grandparent );
        node_proposal_two->prepareProposal();
    }
    if ( &new_grandparent->getChild(0) != grandparent && &new_grandparent->getChild(0) != &grandparent->getChild(0) && &new_grandparent->getChild(0) != &grandparent->getChild(1) )
    {
        left_proposal_two->assignNode( &new_grandparent->getChild(0) );
        left_proposal_two->prepareProposal();
    }
    if ( &new_grandparent->getChild(0) != grandparent && &new_grandparent->getChild(1) != &grandparent->getChild(0) && &new_grandparent->getChild(1) != &grandparent->getChild(1) )
    {
        right_proposal_two->assignNode( &new_grandparent->getChild(1) );
        right_proposal_two->prepareProposal();
    }

    if ( &parent->getChild(0) != grandparent && &parent->getChild(0) != &grandparent->getChild(0) && &parent->getChild(0) != &grandparent->getChild(1) )
    {
        left_proposal_three->assignNode( &parent->getChild(0) );
        left_proposal_three->prepareProposal();
    }
    if ( &parent->getChild(0) != grandparent && &parent->getChild(1) != &grandparent->getChild(0) && &parent->getChild(1) != &grandparent->getChild(1) )
    {
        right_proposal_three->assignNode( &parent->getChild(1) );
        right_proposal_three->prepareProposal();
    }

//    node_proposal_three->assignNode( parent );
//    node_proposal_three->prepareProposal();
    

    // step 4: sample new node states for the affected nodes
//    std::cout << "Storing node states." << std::endl;
    sampleNodeCharacters( grandparent );        // grandparent is the parent of stored_brother
    sampleNodeCharacters( new_grandparent );    // new_grandparent is the parent of parent
    sampleNodeCharacters( parent );             // parent is the parent of stored_new_brother

//    const std::vector<CharacterEvent*>& node_states_one_final = p->getHistory( *grandparent ).getChildCharacters();
//    std::cout << "Proposed states one: ";
//    for (size_t site_index = 0; site_index < num_sites; ++site_index)
//    {
//        size_t s = node_states_one_final[site_index]->getState();
//        std::cout << s << " ";
//    }
//    std::cout << std::endl;
//
//    const std::vector<CharacterEvent*>& node_states_two_final = p->getHistory( *new_grandparent ).getChildCharacters();
//    std::cout << "Proposed states two: ";
//    for (size_t site_index = 0; site_index < num_sites; ++site_index)
//    {
//        size_t s = node_states_two_final[site_index]->getState();
//        std::cout << s << " ";
//    }
//    std::cout << std::endl;
//
//    const std::vector<CharacterEvent*>& node_states_three_final = p->getHistory( *parent ).getChildCharacters();
//    std::cout << "Proposed states three: ";
//    for (size_t site_index = 0; site_index < num_sites; ++site_index)
//    {
//        size_t s = node_states_three_final[site_index]->getState();
//        std::cout << s << " ";
//    }
//    std::cout << std::endl;

    // step 5: sample new branch histories
//    std::cout << "Proposing paths." << std::endl;
    double ln_proposal_probability = 0.0;
    ln_proposal_probability += node_proposal_one->doProposal();
    ln_proposal_probability += left_proposal_one->doProposal();
    ln_proposal_probability += right_proposal_one->doProposal();
    
    if ( &new_grandparent != &grandparent && new_grandparent != &grandparent->getChild(0) && new_grandparent != &grandparent->getChild(1) )
    {
        ln_proposal_probability += node_proposal_two->doProposal();
    }
    if ( &new_grandparent->getChild(0) != grandparent && &new_grandparent->getChild(0) != &grandparent->getChild(0) && &new_grandparent->getChild(0) != &grandparent->getChild(1) )
    {
        ln_proposal_probability += left_proposal_two->doProposal();
    }
    if ( &new_grandparent->getChild(0) != grandparent && &new_grandparent->getChild(1) != &grandparent->getChild(0) && &new_grandparent->getChild(1) != &grandparent->getChild(1) )
    {
        ln_proposal_probability += right_proposal_two->doProposal();
    }
    
    if ( &parent->getChild(0) != grandparent && &parent->getChild(0) != &grandparent->getChild(0) && &parent->getChild(0) != &grandparent->getChild(1) )
    {
        ln_proposal_probability += left_proposal_three->doProposal();
    }
    if ( &parent->getChild(0) != grandparent && &parent->getChild(1) != &grandparent->getChild(0) && &parent->getChild(1) != &grandparent->getChild(1) )
    {
        ln_proposal_probability += right_proposal_three->doProposal();
    }
    
    std::cout << "Proposal ratio: " << ln_proposal_probability << std::endl;
//    std::cout << "Proposed probability: " << p->computeLnProbability() << std::endl;
//    std::cout << "Proposal done." << std::endl;
    
    // step 7: return the hastings ratio
    return 0.0;
    
//    return ln_proposal_probability;
    
}


/**
 *
 */
template<class charType>
void RevBayesCore::FixedNodeheightPruneAndRegraftCharacterHistoryProposal<charType>::prepareProposal( void )
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
template<class charType>
void RevBayesCore::FixedNodeheightPruneAndRegraftCharacterHistoryProposal<charType>::printParameterSummary(std::ostream &o, bool name_only) const
{
    
    // no parameters
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
template<class charType>
void RevBayesCore::FixedNodeheightPruneAndRegraftCharacterHistoryProposal<charType>::undoProposal( void )
{
    
//    std::cout << "Rejecting proposal." << std::endl;    
    
    // we undo the proposal only if it didn't fail
    if ( !failed )
    {
        
//        std::cout << "Rejecting successful proposal." << std::endl;
        
        // now undo changes to the histories
        TreeHistoryCtmc<charType>* p = dynamic_cast< TreeHistoryCtmc<charType>* >( &ctmc->getDistribution() );
        if ( p == NULL )
        {
            throw RbException("Failed cast.");
        }
        size_t num_sites = p->getNumberOfSites();
        const std::vector<BranchHistory*>& histories = p->getHistories();
        
//        std::cout << "Probability before rejection: " << p->computeLnProbability() << std::endl;
        
        // undo the proposal
        TopologyNode& parent = stored_new_brother->getParent();
        TopologyNode& new_grandparent = parent.getParent();
        TopologyNode& grandparent = stored_brother->getParent();
        
        // prune
        new_grandparent.removeChild( &parent );
        parent.removeChild( stored_new_brother );
        new_grandparent.addChild( stored_new_brother );
        stored_new_brother->setParent( &new_grandparent );
        
        // regraft
        grandparent.removeChild( stored_brother );
        parent.addChild( stored_brother );
        stored_brother->setParent( &parent );
        grandparent.addChild( &parent );
        parent.setParent( &grandparent );

        // restore the grandparent states
        std::vector<CharacterEvent*>& node_states_one  = histories[grandparent.getIndex()]->getChildCharacters();
        std::vector<CharacterEvent*>& left_states_one  = histories[grandparent.getChild(0).getIndex()]->getParentCharacters();
        std::vector<CharacterEvent*>& right_states_one = histories[grandparent.getChild(1).getIndex()]->getParentCharacters();
//        std::cout << "Proposed states one: ";
//        for (size_t site_index = 0; site_index < num_sites; ++site_index)
//        {
//            size_t s = node_states_one[site_index]->getState();
//            std::cout << s << " ";
//        }
//        std::cout << std::endl;
        
        // restore the new_grandparent states
        std::vector<CharacterEvent*>& node_states_two  = histories[new_grandparent.getIndex()]->getChildCharacters();
        std::vector<CharacterEvent*>& left_states_two  = histories[new_grandparent.getChild(0).getIndex()]->getParentCharacters();
        std::vector<CharacterEvent*>& right_states_two = histories[new_grandparent.getChild(1).getIndex()]->getParentCharacters();
//        std::cout << "Proposed states two: ";
//        for (size_t site_index = 0; site_index < num_sites; ++site_index)
//        {
//            size_t s = node_states_two[site_index]->getState();
//            std::cout << s << " ";
//        }
//        std::cout << std::endl;
        
        // restore the parent states
        std::vector<CharacterEvent*>& node_states_three  = histories[parent.getIndex()]->getChildCharacters();
        std::vector<CharacterEvent*>& left_states_three  = histories[parent.getChild(0).getIndex()]->getParentCharacters();
        std::vector<CharacterEvent*>& right_states_three = histories[parent.getChild(1).getIndex()]->getParentCharacters();
//        std::cout << "Proposed states three: ";
//        for (size_t site_index = 0; site_index < num_sites; ++site_index)
//        {
//            size_t s = node_states_three[site_index]->getState();
//            std::cout << s << " ";
//        }
//        std::cout << std::endl;
        
//        std::cout << "Restored states one: ";
        for (size_t site_index = 0; site_index < num_sites; ++site_index)
        {
            size_t s = stored_node_states_one[site_index];
            static_cast<CharacterEventDiscrete*>(node_states_one[site_index])->setState(s);
            static_cast<CharacterEventDiscrete*>(left_states_one[site_index])->setState(s);
            static_cast<CharacterEventDiscrete*>(right_states_one[site_index])->setState(s);
//            std::cout << s << " ";
        }
//        std::cout << std::endl;

//        std::cout << "Restored states two: ";
        for (size_t site_index = 0; site_index < num_sites; ++site_index)
        {
            size_t s = stored_node_states_two[site_index];
            static_cast<CharacterEventDiscrete*>(node_states_two[site_index])->setState(s);
            static_cast<CharacterEventDiscrete*>(left_states_two[site_index])->setState(s);
            static_cast<CharacterEventDiscrete*>(right_states_two[site_index])->setState(s);
//            std::cout << s << " ";
        }
//        std::cout << std::endl;

//        std::cout << "Restored states three: ";
        for (size_t site_index = 0; site_index < num_sites; ++site_index)
        {
            size_t s = stored_node_states_three[site_index];
            static_cast<CharacterEventDiscrete*>(node_states_three[site_index])->setState(s);
            static_cast<CharacterEventDiscrete*>(left_states_three[site_index])->setState(s);
            static_cast<CharacterEventDiscrete*>(right_states_three[site_index])->setState(s);
//            std::cout << s << " ";
        }
//        std::cout << std::endl;
        
        // restore path states
//        std::cout << "Restoring branch histories." << std::endl;
        
//        std::cout << "Restoring grandparent histories." << std::endl;
        node_proposal_one->undoProposal();
        left_proposal_one->undoProposal();
        right_proposal_one->undoProposal();

//        std::cout << "Restoring new_grandparent histories." << std::endl;
        if ( &new_grandparent != &grandparent && &new_grandparent != &grandparent.getChild(0) && &new_grandparent != &grandparent.getChild(1) )
        {
            node_proposal_two->undoProposal();
        }
        if ( &new_grandparent.getChild(0) != &grandparent && &new_grandparent.getChild(0) != &grandparent.getChild(0) && &new_grandparent.getChild(0) != &grandparent.getChild(1) )
        {
            left_proposal_two->undoProposal();
        }
        if ( &new_grandparent.getChild(0) != &grandparent && &new_grandparent.getChild(1) != &grandparent.getChild(0) && &new_grandparent.getChild(1) != &grandparent.getChild(1) )
        {
            right_proposal_two->undoProposal();
        }
        
//        std::cout << "Restoring parent histories." << std::endl;
        if ( &parent.getChild(0) != &grandparent && &parent.getChild(0) != &grandparent.getChild(0) && &parent.getChild(0) != &grandparent.getChild(1) )
        {
            left_proposal_three->undoProposal();
        }
        if ( &parent.getChild(0) != &grandparent && &parent.getChild(1) != &grandparent.getChild(0) && &parent.getChild(1) != &grandparent.getChild(1) )
        {
            right_proposal_three->undoProposal();
        }

//        std::cout << "Probability after rejection: " << p->computeLnProbability() << std::endl;
   

    }
    else
    {
//        std::cout << "Rejecting failed proposal." << std::endl;
    }
    
}

template<class charType>
void RevBayesCore::FixedNodeheightPruneAndRegraftCharacterHistoryProposal<charType>::sampleNodeCharacters( TopologyNode* node )
{
    
    if ( node->isTip() == false )
    {
        TreeHistoryCtmc<charType>* c = dynamic_cast< TreeHistoryCtmc<charType>* >(&ctmc->getDistribution());
        if ( ctmc == NULL )
        {
            throw RbException("Failed cast.");
        }
        size_t num_sites = c->getNumberOfSites();
        const std::vector<BranchHistory*>& histories = c->getHistories();
        
        TopologyNode &left_child  = node->getChild(0);
        TopologyNode &right_child = node->getChild(1);
        
        double node_age   = node->getAge();
        double left_age   = left_child.getAge();
        double right_age  = right_child.getAge();
        
        double node_rate  = c->getBranchRate( node->getIndex()       );
        double left_rate  = c->getBranchRate( left_child.getIndex()  );
        double right_rate = c->getBranchRate( right_child.getIndex() );
        
        // get transition probs
        const RateGenerator& rm = rate_generator->getValue();
        rm.calculateTransitionProbabilities(node_age, left_age,  left_rate,  left_tp_matrix);
        rm.calculateTransitionProbabilities(node_age, right_age, right_rate, right_tp_matrix);
        
        // states for conditional sampling probs
        const std::vector<CharacterEvent*>& leftChildState  = histories[left_child.getIndex()]->getChildCharacters();
        const std::vector<CharacterEvent*>& rightChildState = histories[right_child.getIndex()]->getChildCharacters();
        
        // states to update
        std::vector<CharacterEvent*>& nodeChildState   = histories[node->getIndex()]->getChildCharacters();
        std::vector<CharacterEvent*>& leftParentState  = histories[left_child.getIndex()]->getParentCharacters();
        std::vector<CharacterEvent*>& rightParentState = histories[right_child.getIndex()]->getParentCharacters();
        
        if ( node->isRoot() == false )
        {
            
            double parent_age = node->getParent().getAge();
            rm.calculateTransitionProbabilities(parent_age, node_age, node_rate, node_tp_matrix);
            
            const std::vector<CharacterEvent*>& nodeParentState = histories[node->getIndex()]->getParentCharacters();
            
            for (size_t site_index = 0; site_index < num_sites; ++site_index)
            {
                size_t ancS  = static_cast<CharacterEventDiscrete*>(nodeParentState[site_index])->getState();
                size_t desS1 = static_cast<CharacterEventDiscrete*>(leftChildState[site_index])->getState();
                size_t desS2 = static_cast<CharacterEventDiscrete*>(rightChildState[site_index])->getState();
                
                double u = GLOBAL_RNG->uniform01();
                
                std::vector<double> state_probs(num_states, 0.0);
                double sum = 0.0;
                for ( size_t i=0; i < num_states; ++i )
                {
                    double p = node_tp_matrix[ancS][i] * left_tp_matrix[i][desS1] * right_tp_matrix[i][desS2];
                    sum += p;
                    state_probs[i] = p;
                }
                
                unsigned int s = 0;
                for ( size_t i=0; i < num_states; ++i )
                {
                    u -= (state_probs[i]/sum);
                    if (u <= 0.0 )
                    {
                        //                        s = i;
                        break;
                    }
                    ++s;
                }
                
                static_cast<CharacterEventDiscrete*>(nodeChildState[site_index])->setState(s);
                static_cast<CharacterEventDiscrete*>(leftParentState[site_index])->setState(s);
                static_cast<CharacterEventDiscrete*>(rightParentState[site_index])->setState(s);
            }
            
        }
        else
        {

            std::vector<double> rf = c->getRootFrequencies();
            
            for (size_t site_index = 0; site_index < num_sites; ++site_index)
            {
                size_t desS1 = static_cast<CharacterEventDiscrete*>(leftChildState[site_index])->getState();
                size_t desS2 = static_cast<CharacterEventDiscrete*>(rightChildState[site_index])->getState();
                
                double u = GLOBAL_RNG->uniform01();
                
                std::vector<double> state_probs(num_states, 0.0);
                double sum = 0.0;
                for ( size_t i=0; i < num_states; ++i )
                {
                    double p = rf[i] * left_tp_matrix[i][desS1] * right_tp_matrix[i][desS2];
                    sum += p;
                    state_probs[i] = p;
                }
                
                unsigned int s = 0;
                for ( size_t i=0; i < num_states; ++i )
                {
                    u -= (state_probs[i]/sum);
                    if (u <= 0.0 )
                    {
                        //                        s = i;
                        break;
                    }
                    ++s;
                }
                
                static_cast<CharacterEventDiscrete*>(nodeChildState[site_index])->setState(s);
                static_cast<CharacterEventDiscrete*>(leftParentState[site_index])->setState(s);
                static_cast<CharacterEventDiscrete*>(rightParentState[site_index])->setState(s);
                
            }
            
        }
        
    }
    
}

template<class charType>
void RevBayesCore::FixedNodeheightPruneAndRegraftCharacterHistoryProposal<charType>::setRateGenerator(const TypedDagNode<RateGenerator> *d)
{
    
    rate_generator = d;
    num_states = rate_generator->getValue().getNumberOfStates();
    
    node_proposal_one->setRateGenerator( rate_generator );
    left_proposal_one->setRateGenerator( rate_generator );
    right_proposal_one->setRateGenerator( rate_generator );
    node_proposal_two->setRateGenerator( rate_generator );
    left_proposal_two->setRateGenerator( rate_generator );
    right_proposal_two->setRateGenerator( rate_generator );
//    node_proposal_three->setRateGenerator( rate_generator );
    left_proposal_three->setRateGenerator( rate_generator );
    right_proposal_three->setRateGenerator( rate_generator );
    
    node_tp_matrix  = TransitionProbabilityMatrix( num_states );
    left_tp_matrix  = TransitionProbabilityMatrix( num_states );
    right_tp_matrix = TransitionProbabilityMatrix( num_states );
    
}

template<class charType>
void RevBayesCore::FixedNodeheightPruneAndRegraftCharacterHistoryProposal<charType>::setRateGenerator(const TypedDagNode<RateGeneratorSequence> *d)
{
    
    rate_generator_sequence = d;
    num_states = rate_generator_sequence->getValue().getNumberOfStates();
    
    node_proposal_one->setRateGenerator( rate_generator_sequence );
    left_proposal_one->setRateGenerator( rate_generator_sequence );
    right_proposal_one->setRateGenerator( rate_generator_sequence );
    node_proposal_two->setRateGenerator( rate_generator_sequence );
    left_proposal_two->setRateGenerator( rate_generator_sequence );
    right_proposal_two->setRateGenerator( rate_generator_sequence );
//    node_proposal_three->setRateGenerator( rate_generator_sequence );
    left_proposal_three->setRateGenerator( rate_generator_sequence );
    right_proposal_three->setRateGenerator( rate_generator_sequence );
    
    node_tp_matrix  = TransitionProbabilityMatrix( num_states );
    left_tp_matrix  = TransitionProbabilityMatrix( num_states );
    right_tp_matrix = TransitionProbabilityMatrix( num_states );
    
}

/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
template<class charType>
void RevBayesCore::FixedNodeheightPruneAndRegraftCharacterHistoryProposal<charType>::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    if ( oldN == tree )
    {
        tree = static_cast<StochasticNode<Tree>* >(newN);
    }
    else if ( oldN == ctmc )
    {
        ctmc = static_cast<StochasticNode<AbstractHomologousDiscreteCharacterData>* >(newN);
    }
    else if ( oldN == rate_generator )
    {
        rate_generator = static_cast<TypedDagNode<RateGenerator>* >(newN);
    }
    else if ( oldN == rate_generator_sequence )
    {
        rate_generator_sequence = static_cast<DeterministicNode<RateGeneratorSequence>* >(newN);
    }
    
    node_proposal_one->swapNodeInternal(oldN, newN);
    left_proposal_one->swapNodeInternal(oldN, newN);
    right_proposal_one->swapNodeInternal(oldN, newN);
    
    node_proposal_two->swapNodeInternal(oldN, newN);
    left_proposal_two->swapNodeInternal(oldN, newN);
    right_proposal_two->swapNodeInternal(oldN, newN);
   
//    node_proposal_three->swapNodeInternal(oldN, newN);
    left_proposal_three->swapNodeInternal(oldN, newN);
    right_proposal_three->swapNodeInternal(oldN, newN);
    
}

template<class charType>
void RevBayesCore::FixedNodeheightPruneAndRegraftCharacterHistoryProposal<charType>::setProposalTuningParameter(double tp)
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
template<class charType>
void RevBayesCore::FixedNodeheightPruneAndRegraftCharacterHistoryProposal<charType>::tune( double rate )
{
    
    // nothing to tune
    
}


#endif

