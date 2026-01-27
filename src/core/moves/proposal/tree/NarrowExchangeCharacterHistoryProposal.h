#ifndef NarrowExchangeCharacterHistoryProposal_H
#define NarrowExchangeCharacterHistoryProposal_H

#include "BranchHistoryDiscrete.h"
#include "CharacterEventDiscrete.h"
#include "HomologousDiscreteCharacterData.h"
#include "PathRejectionSampleProposal.h"
#include "NodeRejectionSampleProposal.h"
#include "Proposal.h"
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
     * The narrow-exchange operator.
     *
     * A narrow-exchange proposal is a NNI (nearest neighbour interchange) proposal on rooted trees without changing the node age.
     * That is, we pick a random node which is not the root and neither its parent is the root.
     * Then, we try to exchange the picked node with its uncle. This move will automatically fail if the uncle is older than the parent.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2012-07-12, version 1.0
     *
     */
    
    template<class charType>
    class NarrowExchangeCharacterHistoryProposal : public Proposal {

    public:
        NarrowExchangeCharacterHistoryProposal( StochasticNode<Tree> *t, StochasticNode<AbstractHomologousDiscreteCharacterData>* c);                                               //!<  constructor

        // Basic utility functions
        void                                                        cleanProposal(void);                                            //!< Clean up proposal
        NarrowExchangeCharacterHistoryProposal*                     clone(void) const;                                              //!< Clone object
        double                                                      computeLnProposal( TopologyNode* node );
        double                                                      computeNodeStateProbability( TopologyNode* parent);
        double                                                      doProposal(void);                                               //!< Perform proposal
        const std::string&                                          getProposalName(void) const;                                    //!< Get the name of the proposal for summary printing
        double                                                      getProposalTuningParameter(void) const;
        void                                                        prepareProposal(void);                                          //!< Prepare the proposal
        void                                                        printParameterSummary(std::ostream &o, bool name_only) const;                   //!< Print the parameter summary
        void                                                        sampleNodeCharacters(TopologyNode* node);                       //!< Sample the characters at the node
        void                                                        sampleNodeCharactersJoint(TopologyNode* parent);                                         //!< Sample the characters at the node
        void                                                        setRateGenerator(const TypedDagNode<RateGenerator> *d);
        void                                                        setRateGenerator(const TypedDagNode<RateGeneratorSequence> *d);
        void                                                        setProposalTuningParameter(double tp);
        void                                                        tune(double r);                                                 //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                                        undoProposal(void);                                             //!< Reject the proposal

    protected:

        void                                                        swapNodeInternal(DagNode *oldN, DagNode *newN);                 //!< Swap the DAG nodes on which the Proposal is working on

    private:

        // parameters
        StochasticNode<Tree>*                                       tree;                                                           //!< The tree the Proposal is working on
        StochasticNode<AbstractHomologousDiscreteCharacterData>*    ctmc;                                                           //!< The character history the Proposal is working on
        const TypedDagNode<RateGenerator>*                          rate_generator;
        const TypedDagNode<RateGeneratorSequence>*                  rate_generator_sequence;

        // dimensions
        size_t                                                      num_states;
        
        // proposals
        PathRejectionSampleProposal<charType>*                      node_proposal;
        PathRejectionSampleProposal<charType>*                      parent_proposal;
        PathRejectionSampleProposal<charType>*                      brother_proposal;
        PathRejectionSampleProposal<charType>*                      uncle_proposal;
        PathRejectionSampleProposal<charType>*                      grandparent_proposal;
        
        // transition probability matrices
        TransitionProbabilityMatrix                                 node_tp_matrix;
        TransitionProbabilityMatrix                                 brother_tp_matrix;
        TransitionProbabilityMatrix                                 uncle_tp_matrix;
        TransitionProbabilityMatrix                                 parent_tp_matrix;
        TransitionProbabilityMatrix                                 grandparent_tp_matrix;
        
        // stored objects to undo proposal
        bool                                                        failed;
        TopologyNode*                                               stored_chosen_node;
        TopologyNode*                                               stored_uncle;
        
        std::vector<size_t>                                         stored_parent_node_states;
        std::vector<size_t>                                         stored_grandparent_node_states;
        
        size_t                                                      stored_index_parent;
        size_t                                                      stored_index_node;
        size_t                                                      stored_index_uncle;
        size_t                                                      stored_index_brother;
        size_t                                                      proposed_index_parent;
        size_t                                                      proposed_index_node;
        size_t                                                      proposed_index_uncle;
        size_t                                                      proposed_index_brother;
        
    };

}

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
template<class charType>
RevBayesCore::NarrowExchangeCharacterHistoryProposal<charType>::NarrowExchangeCharacterHistoryProposal( StochasticNode<Tree> *t, StochasticNode<AbstractHomologousDiscreteCharacterData>* c ) : Proposal(),
    tree( t ),
    ctmc( c ),
    rate_generator( NULL ),
    rate_generator_sequence( NULL ),
    num_states( 2 ),
    node_tp_matrix( 2 ),
    brother_tp_matrix( 2 ),
    uncle_tp_matrix( 2 ),
    parent_tp_matrix( 2 ),
    grandparent_tp_matrix( 2 )
{
    
    // tell the base class to add the nodes
    addNode( tree );
    addNode( ctmc );
    
    // make the proposals
    node_proposal        = new PathRejectionSampleProposal<charType>( c );
    parent_proposal      = new PathRejectionSampleProposal<charType>( c );
    brother_proposal     = new PathRejectionSampleProposal<charType>( c );
    uncle_proposal       = new PathRejectionSampleProposal<charType>( c );
    grandparent_proposal = new PathRejectionSampleProposal<charType>( c );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
template<class charType>
void RevBayesCore::NarrowExchangeCharacterHistoryProposal<charType>::cleanProposal( void )
{
    node_proposal->cleanProposal();
    parent_proposal->cleanProposal();
    brother_proposal->cleanProposal();
    uncle_proposal->cleanProposal();
    grandparent_proposal->cleanProposal();
    
    stored_parent_node_states.clear();
    stored_grandparent_node_states.clear();
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
template<class charType>
RevBayesCore::NarrowExchangeCharacterHistoryProposal<charType>* RevBayesCore::NarrowExchangeCharacterHistoryProposal<charType>::clone( void ) const
{
    
    return new NarrowExchangeCharacterHistoryProposal( *this );
}

template<class charType>
double RevBayesCore::NarrowExchangeCharacterHistoryProposal<charType>::computeNodeStateProbability( TopologyNode* parent )
{
    double lnP = 0.0;
    
    // retrieve the character history object
    TreeHistoryCtmc<charType>* c = dynamic_cast< TreeHistoryCtmc<charType>* >(&ctmc->getDistribution());
    if ( ctmc == NULL )
    {
        throw RbException("Failed cast.");
    }
    size_t num_sites = c->getNumberOfSites();
    const std::vector<BranchHistory*>& histories = c->getHistories();
    
    // get all the nodes
    TopologyNode& grandparent = parent->getParent();
    TopologyNode* uncle = &grandparent.getChild( 0 );
    // check if we got the correct child
    if ( uncle == parent )
    {
        uncle = &grandparent.getChild( 1 );
    }
    TopologyNode* node    = &parent->getChild( 0 );
    TopologyNode* brother = &parent->getChild( 1 );
    
    // states for conditional sampling probs
    const std::vector<CharacterEvent*>& nodeChildState    = histories[node->getIndex()]->getChildCharacters();
    const std::vector<CharacterEvent*>& brotherChildState = histories[brother->getIndex()]->getChildCharacters();
    const std::vector<CharacterEvent*>& uncleChildState   = histories[uncle->getIndex()]->getChildCharacters();
    
    const std::vector<CharacterEvent*>& parentChildState      = histories[parent->getIndex()]->getChildCharacters();
    const std::vector<CharacterEvent*>& grandparentChildState = histories[grandparent.getIndex()]->getChildCharacters();
    
    // update the transition probs
    const RateGenerator& rm = rate_generator->getValue();
    rm.calculateTransitionProbabilities(node->getParent().getAge(),    node->getAge(),    c->getBranchRate( node->getIndex() ),    node_tp_matrix);
    rm.calculateTransitionProbabilities(brother->getParent().getAge(), brother->getAge(), c->getBranchRate( brother->getIndex() ), brother_tp_matrix);
    rm.calculateTransitionProbabilities(uncle->getParent().getAge(),   uncle->getAge(),   c->getBranchRate( uncle->getIndex() ),   uncle_tp_matrix);
    rm.calculateTransitionProbabilities(parent->getParent().getAge(),  parent->getAge(),  c->getBranchRate( parent->getIndex() ),  parent_tp_matrix);
    
    if ( grandparent.isRoot() == false )
    {
        
        // update the transition probs
        rm.calculateTransitionProbabilities(grandparent.getParent().getAge(), grandparent.getAge(),  c->getBranchRate( grandparent.getIndex() ), grandparent_tp_matrix);
        
        // get the state of the great-grandparent
        const std::vector<CharacterEvent*>& greatGrandparentState = histories[grandparent.getIndex()]->getParentCharacters();
        
        // for each site...
        for (size_t site_index = 0; site_index < num_sites; ++site_index)
        {
            
            size_t node_state              = static_cast<CharacterEventDiscrete*>(nodeChildState[site_index])->getState();
            size_t brother_state           = static_cast<CharacterEventDiscrete*>(brotherChildState[site_index])->getState();
            size_t uncle_state             = static_cast<CharacterEventDiscrete*>(uncleChildState[site_index])->getState();
            size_t great_grandparent_state = static_cast<CharacterEventDiscrete*>(greatGrandparentState[site_index])->getState();
            
            double sum = 0.0;
            
            std::vector<double> state_probs(num_states * num_states);
            
            // for each grandparent state...
            for (size_t grandparent_state = 0; grandparent_state < num_states; ++grandparent_state)
            {
                
                // for each parent state...
                for (size_t parent_state = 0; parent_state < num_states; ++parent_state)
                {
                    double p = parent_tp_matrix[grandparent_state][parent_state] * node_tp_matrix[parent_state][node_state] * brother_tp_matrix[parent_state][brother_state] * uncle_tp_matrix[grandparent_state][uncle_state] * grandparent_tp_matrix[great_grandparent_state][grandparent_state];
                    sum += p;
                    state_probs[grandparent_state * num_states + parent_state] = p;
                }
                
            }
            
            size_t parent_state = static_cast<CharacterEventDiscrete*>(parentChildState[site_index])->getState();
            size_t grandparent_state = static_cast<CharacterEventDiscrete*>(grandparentChildState[site_index])->getState();
            lnP += log(state_probs[ num_states * grandparent_state + parent_state ] / sum);

        }
        
    }
    else
    {
        
        // update the transition probs
        std::vector<double> rf = c->getRootFrequencies();
        
        // for each site...
        for (size_t site_index = 0; site_index < num_sites; ++site_index)
        {
            
            size_t node_state              = static_cast<CharacterEventDiscrete*>(nodeChildState[site_index])->getState();
            size_t brother_state           = static_cast<CharacterEventDiscrete*>(brotherChildState[site_index])->getState();
            size_t uncle_state             = static_cast<CharacterEventDiscrete*>(uncleChildState[site_index])->getState();
            
            double sum = 0.0;
            
            std::vector<double> state_probs(num_states * num_states);
            
            // for each grandparent state...
            for (size_t grandparent_state = 0; grandparent_state < num_states; ++grandparent_state)
            {
                
                // for each parent state...
                for (size_t parent_state = 0; parent_state < num_states; ++parent_state)
                {
                    double p = parent_tp_matrix[grandparent_state][parent_state] * node_tp_matrix[parent_state][node_state] * brother_tp_matrix[parent_state][brother_state] * uncle_tp_matrix[grandparent_state][uncle_state] * rf[grandparent_state];
                    sum += p;
                    state_probs[grandparent_state * num_states + parent_state] = p;
                }
                
            }
            
            size_t parent_state = static_cast<CharacterEventDiscrete*>(parentChildState[site_index])->getState();
            size_t grandparent_state = static_cast<CharacterEventDiscrete*>(grandparentChildState[site_index])->getState();
            lnP += log(state_probs[ num_states * grandparent_state + parent_state ] / sum);

        }
        
    }

    
    return lnP;
}

template<class charType>
double RevBayesCore::NarrowExchangeCharacterHistoryProposal<charType>::computeLnProposal( TopologyNode* node )
{
    double lnP = 0.0;
    
    if ( node->isTip() == false )
    {
        TreeHistoryCtmc<charType>* c = dynamic_cast< TreeHistoryCtmc<charType>* >(&ctmc->getDistribution());
        if ( ctmc == NULL )
        {
            throw RbException("Failed cast.");
        }
        
        size_t num_sites = c->getNumberOfSites();
        
        std::vector<BranchHistory*> histories = c->getHistories();
        
        TopologyNode &left_child  = node->getChild(0);
        TopologyNode &right_child = node->getChild(1);
        
        double node_age     = node->getAge();
        double left_age     = left_child.getAge();
        double right_age    = right_child.getAge();
        
        double node_rate  = c->getBranchRate( node->getIndex() );
        double left_rate  = c->getBranchRate( left_child.getIndex() );
        double right_rate = c->getBranchRate( right_child.getIndex() );
        
        // get transition probs
        const RateGenerator& rm = rate_generator->getValue();
        rm.calculateTransitionProbabilities(node_age, left_age,  left_rate,  brother_tp_matrix);
        rm.calculateTransitionProbabilities(node_age, right_age, right_rate, uncle_tp_matrix);
        
        // states for conditional sampling probs
        const std::vector<CharacterEvent*>& leftChildState  = histories[node->getChild(0).getIndex()]->getChildCharacters();
        const std::vector<CharacterEvent*>& rightChildState = histories[node->getChild(1).getIndex()]->getChildCharacters();
        
        // states to update
        const std::vector<CharacterEvent*>& nodeChildState  = histories[node->getIndex()]->getChildCharacters();
        
        if ( node->isRoot() == false )
        {
            double parent_age   = node->getParent().getAge();
            rm.calculateTransitionProbabilities(parent_age, node_age, node_rate, node_tp_matrix);
            
            const std::vector<CharacterEvent*>& nodeParentState = histories[node->getIndex()]->getParentCharacters();
            
            for (size_t site_index = 0; site_index < num_sites; ++site_index)
            {
                size_t ancS  = static_cast<CharacterEventDiscrete*>(nodeParentState[site_index])->getState();
                size_t desS1 = static_cast<CharacterEventDiscrete*>(leftChildState[site_index])->getState();
                size_t desS2 = static_cast<CharacterEventDiscrete*>(rightChildState[site_index])->getState();
                
                std::vector<double> state_probs(num_states,0.0);
                double sum = 0.0;
                for ( size_t i=0; i < num_states; ++i )
                {
                    double p = node_tp_matrix[ancS][i] * brother_tp_matrix[i][desS1] * uncle_tp_matrix[i][desS2];
                    sum += p;
                    state_probs[i] = p;
                }
                
                size_t s = static_cast<CharacterEventDiscrete*>(nodeChildState[site_index])->getState();
                lnP += log(state_probs[s]/sum);
            }
            
        }
        else
        {
            
            std::vector<double> rf = c->getRootFrequencies();
            
            for (size_t site_index = 0; site_index < num_sites; ++site_index)
            {
                size_t desS1 = static_cast<CharacterEventDiscrete*>(leftChildState[site_index])->getState();
                size_t desS2 = static_cast<CharacterEventDiscrete*>(rightChildState[site_index])->getState();
                
                std::vector<double> state_probs(num_states,0.0);
                double sum = 0.0;
                for ( size_t i=0; i < num_states; ++i )
                {
                    double p = rf[i] * brother_tp_matrix[i][desS1] * uncle_tp_matrix[i][desS2];
                    sum += p;
                    state_probs[i] = p;
                }
                
                size_t s = static_cast<CharacterEventDiscrete*>(nodeChildState[site_index])->getState();
                lnP += log(state_probs[s]/sum);
                
            }
            
        }
        
    }
    
    return lnP;
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
template<class charType>
const std::string& RevBayesCore::NarrowExchangeCharacterHistoryProposal<charType>::getProposalName( void ) const
{
    static std::string name = "NarrowExchangeCharacterHistoryProposal";
    
    return name;
}


template<class charType>
double RevBayesCore::NarrowExchangeCharacterHistoryProposal<charType>::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}


/**
 * Perform the proposal.
 *
 * A Beta-simplex proposal randomly changes some values of a simplex, although the other values
 * change too because of the renormalization.
 * First, some random indices are drawn. Then, the proposal draws a new simplex
 *   u ~ Beta(val[index] * alpha)
 * where alpha is the tuning parameter.The new value is set to u.
 * The simplex is then renormalized.
 *
 * \return The hastings ratio.
 */
template<class charType>
double RevBayesCore::NarrowExchangeCharacterHistoryProposal<charType>::doProposal( void )
{
    
    // Get random number generator
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    Tree& tau = tree->getValue();
    
    // pick a random node which is not the root and neither a direct descendant of the root
    TopologyNode* node;
    do {
        double u = rng->uniform01();
        size_t index = size_t( std::floor(tau.getNumberOfNodes() * u) );
        node = &tau.getNode(index);
    } while ( node->isRoot() || node->getParent().isRoot() );
    
    TopologyNode& parent      = node->getParent();
    TopologyNode& grandparent = parent.getParent();
    TopologyNode* uncle       = &grandparent.getChild( 0 );
    // check if we got the correct child
    if ( uncle == &parent )
    {
        uncle = &grandparent.getChild( 1 );
    }
    TopologyNode* brother = &parent.getChild( 0 );
    // check if we got the correct child
    if ( brother == node )
    {
        brother = &parent.getChild( 1 );
    }
    
    // we need to work with the times
    double parent_age = parent.getAge();
    double uncle_age  = uncle->getAge();
    
    if ( uncle_age < parent_age )
    {
        failed = false;
        
        // step 0: retrieve the character history
        TreeHistoryCtmc<charType>* p = dynamic_cast< TreeHistoryCtmc<charType>* >( &ctmc->getDistribution() );
        if ( p == NULL )
        {
            throw RbException("Failed cast.");
        }
        size_t num_sites = p->getNumberOfSites();

//        double before = p->computeLnProbability();
        
        // step 1: store all the necessary values
        stored_chosen_node = node;
        stored_uncle       = uncle;

        stored_parent_node_states.clear();
        const std::vector<CharacterEvent*>& parent_node_states = p->getHistory( parent ).getChildCharacters();
        stored_parent_node_states.resize(num_sites, 0);
        for (size_t site_index = 0; site_index < num_sites; ++site_index)
        {
            size_t s = static_cast<CharacterEventDiscrete*>(parent_node_states[site_index])->getState();
            stored_parent_node_states[site_index] = s;
        }
        
        stored_grandparent_node_states.clear();
        const std::vector<CharacterEvent*>& grandparent_node_states = p->getHistory( grandparent ).getChildCharacters();
        stored_grandparent_node_states.resize(num_sites, 0);
        for (size_t site_index = 0; site_index < num_sites; ++site_index)
        {
            size_t s = static_cast<CharacterEventDiscrete*>(grandparent_node_states[site_index])->getState();
            stored_grandparent_node_states[site_index] = s;
        }

        // step 2: assign the nodes for the path proposals
        node_proposal->assignNode( node );
        parent_proposal->assignNode( &parent );
        brother_proposal->assignNode( brother );
        uncle_proposal->assignNode( uncle );
        grandparent_proposal->assignNode( &grandparent );

        // step 3: prepare the path proposals
        node_proposal->prepareProposal();
        parent_proposal->prepareProposal();
        brother_proposal->prepareProposal();
        uncle_proposal->prepareProposal();
        grandparent_proposal->prepareProposal();

        // step 4: exchange the two nodes
        grandparent.removeChild( uncle );
        parent.removeChild( node );
        grandparent.addChild( node );
        parent.addChild( uncle );
        node->setParent( &grandparent );
        uncle->setParent( &parent );

        // step 5: sample new node states for the parent and grandparent
        sampleNodeCharactersJoint( &parent );

        // step 6: sample new branch histories
        double ln_proposal_probability = 0.0;
        ln_proposal_probability += node_proposal->doProposal();
        ln_proposal_probability += parent_proposal->doProposal();
        ln_proposal_probability += brother_proposal->doProposal();
        ln_proposal_probability += uncle_proposal->doProposal();
        ln_proposal_probability += grandparent_proposal->doProposal();
        
//        double after = p->computeLnProbability();
        
        // step 7: return the hastings ratio
        return 0.0;
//        return ln_proposal_probability - before + after;
        
    }
    else
    {
        failed = true;
        return RbConstants::Double::neginf;
    }
    
}

/**
 *
 */
template<class charType>
void RevBayesCore::NarrowExchangeCharacterHistoryProposal<charType>::prepareProposal( void )
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
void RevBayesCore::NarrowExchangeCharacterHistoryProposal<charType>::printParameterSummary(std::ostream &o, bool name_only) const
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
void RevBayesCore::NarrowExchangeCharacterHistoryProposal<charType>::undoProposal( void )
{
    
//    std::cout << "Undoing proposal." << std::endl;

    // we undo the proposal only if it didn't fail
    if ( failed == false )
    {
        
//        std::cout << "Rejecting proposal." << std::endl;
        
        // undo the proposal
        TopologyNode& parent      = stored_uncle->getParent();
        TopologyNode& grandparent = stored_chosen_node->getParent();
        TopologyNode& uncle       = *stored_uncle;
        TopologyNode* brother     = &parent.getChild( 0 );
        // check if we got the correct child
        if ( brother == &uncle )
        {
            brother = &parent.getChild( 1 );
        }
        
        // now exchange the two nodes
        grandparent.removeChild( stored_chosen_node  );
        parent.removeChild( stored_uncle );
        grandparent.addChild( stored_uncle );
        parent.addChild( stored_chosen_node  );
        stored_uncle->setParent( &grandparent );
        stored_chosen_node->setParent( &parent );
        
        // now undo changes to the histories
        TreeHistoryCtmc<charType>* p = dynamic_cast< TreeHistoryCtmc<charType>* >( &ctmc->getDistribution() );
        if ( p == NULL )
        {
            throw RbException("Failed cast.");
        }
        size_t num_sites = p->getNumberOfSites();
        CharacterHistoryDiscrete& histories = p->getHistories();
        
        // restore the parent states
        std::vector<CharacterEvent*>& parentState     = histories[parent.getIndex()].getChildCharacters();
        std::vector<CharacterEvent*>& leftChildState  = histories[parent.getChild(0).getIndex()].getParentCharacters();
        std::vector<CharacterEvent*>& rightChildState = histories[parent.getChild(1).getIndex()].getParentCharacters();
        for (size_t site_index = 0; site_index < num_sites; ++site_index)
        {
            size_t s = stored_parent_node_states[site_index];
            static_cast<CharacterEventDiscrete*>(parentState[site_index])->setState(s);
            static_cast<CharacterEventDiscrete*>(leftChildState[site_index])->setState(s);
            static_cast<CharacterEventDiscrete*>(rightChildState[site_index])->setState(s);
        }
        
        // restore the grandparent states
        std::vector<CharacterEvent*>& grandparentState           = histories[grandparent.getIndex()].getChildCharacters();
        std::vector<CharacterEvent*>& leftGrandparentChildState  = histories[grandparent.getChild(0).getIndex()].getParentCharacters();
        std::vector<CharacterEvent*>& rightGrandparentChildState = histories[grandparent.getChild(1).getIndex()].getParentCharacters();
        for (size_t site_index = 0; site_index < num_sites; ++site_index)
        {
            size_t s = stored_grandparent_node_states[site_index];
            static_cast<CharacterEventDiscrete*>(grandparentState[site_index])->setState(s);
            static_cast<CharacterEventDiscrete*>(leftGrandparentChildState[site_index])->setState(s);
            static_cast<CharacterEventDiscrete*>(rightGrandparentChildState[site_index])->setState(s);
        }
        
        // restore path states
        node_proposal->undoProposal();
        parent_proposal->undoProposal();
        uncle_proposal->undoProposal();
        brother_proposal->undoProposal();
        grandparent_proposal->undoProposal();
        
    }
    
}

template<class charType>
void RevBayesCore::NarrowExchangeCharacterHistoryProposal<charType>::sampleNodeCharactersJoint( TopologyNode* parent )
{
    
    // retrieve the character history object
    TreeHistoryCtmc<charType>* c = dynamic_cast< TreeHistoryCtmc<charType>* >(&ctmc->getDistribution());
    if ( ctmc == NULL )
    {
        throw RbException("Failed cast.");
    }
    size_t num_sites = c->getNumberOfSites();
    CharacterHistoryDiscrete& histories = c->getHistories();
    
    // get all the nodes
    TopologyNode& grandparent = parent->getParent();
    TopologyNode* uncle = &grandparent.getChild( 0 );
    // check if we got the correct child
    if ( uncle == parent )
    {
        uncle = &grandparent.getChild( 1 );
    }
    TopologyNode* node    = &parent->getChild( 0 );
    TopologyNode* brother = &parent->getChild( 1 );

    // states to update
    std::vector<CharacterEvent*>& nodeParentState        = histories[node->getIndex()].getParentCharacters();
    std::vector<CharacterEvent*>& brotherParentState     = histories[brother->getIndex()].getParentCharacters();
    std::vector<CharacterEvent*>& uncleParentState       = histories[uncle->getIndex()].getParentCharacters();
    std::vector<CharacterEvent*>& parentChildState       = histories[parent->getIndex()].getChildCharacters();
    std::vector<CharacterEvent*>& parentParentState      = histories[parent->getIndex()].getParentCharacters();
    std::vector<CharacterEvent*>& grandparentChildState  = histories[grandparent.getIndex()].getChildCharacters();
    
    // states for conditional sampling probs
    const std::vector<CharacterEvent*>& nodeChildState    = histories[node->getIndex()].getChildCharacters();
    const std::vector<CharacterEvent*>& brotherChildState = histories[brother->getIndex()].getChildCharacters();
    const std::vector<CharacterEvent*>& uncleChildState   = histories[uncle->getIndex()].getChildCharacters();
    
    // update the transition probs
    const RateGenerator& rm = rate_generator->getValue();
    rm.calculateTransitionProbabilities(node->getParent().getAge(),    node->getAge(),    c->getBranchRate( node->getIndex() ),    node_tp_matrix);
    rm.calculateTransitionProbabilities(brother->getParent().getAge(), brother->getAge(), c->getBranchRate( brother->getIndex() ), brother_tp_matrix);
    rm.calculateTransitionProbabilities(uncle->getParent().getAge(),   uncle->getAge(),   c->getBranchRate( uncle->getIndex() ),   uncle_tp_matrix);
    rm.calculateTransitionProbabilities(parent->getParent().getAge(),  parent->getAge(),  c->getBranchRate( parent->getIndex() ),  parent_tp_matrix);
    
    if ( grandparent.isRoot() == false )
    {
        
        // update the transition probs
        rm.calculateTransitionProbabilities(grandparent.getParent().getAge(), grandparent.getAge(),  c->getBranchRate( grandparent.getIndex() ), grandparent_tp_matrix);
        
        // get the state of the great-grandparent
        const std::vector<CharacterEvent*>& greatGrandparentState = histories[grandparent.getIndex()].getParentCharacters();

        // for each site...
        for (size_t site_index = 0; site_index < num_sites; ++site_index)
        {
         
            size_t node_state              = static_cast<CharacterEventDiscrete*>(nodeChildState[site_index])->getState();
            size_t brother_state           = static_cast<CharacterEventDiscrete*>(brotherChildState[site_index])->getState();
            size_t uncle_state             = static_cast<CharacterEventDiscrete*>(uncleChildState[site_index])->getState();
            size_t great_grandparent_state = static_cast<CharacterEventDiscrete*>(greatGrandparentState[site_index])->getState();

            double sum = 0.0;
            
            std::vector<double> state_probs(num_states * num_states);
            
            // for each grandparent state...
            for (size_t grandparent_state = 0; grandparent_state < num_states; ++grandparent_state)
            {
                
                // for each parent state...
                for (size_t parent_state = 0; parent_state < num_states; ++parent_state)
                {
                    double p = parent_tp_matrix[grandparent_state][parent_state] * node_tp_matrix[parent_state][node_state] * brother_tp_matrix[parent_state][brother_state] * uncle_tp_matrix[grandparent_state][uncle_state] * grandparent_tp_matrix[great_grandparent_state][grandparent_state];
                    sum += p;
                    state_probs[grandparent_state * num_states + parent_state] = p;
                }
                
            }
            
            // generate a random uniform number
            double u = GLOBAL_RNG->uniform01();

            size_t new_grandparent_state = 0;
            size_t new_parent_state      = 0;
            
            // sample a new state
            for (size_t grandparent_state = 0; grandparent_state < num_states; ++grandparent_state)
            {
                for (size_t parent_state = 0; parent_state < num_states; ++parent_state)
                {
                    u -= (state_probs[grandparent_state * num_states + parent_state] / sum);
                    if ( u < 0.0 )
                    {
                        new_grandparent_state = grandparent_state;
                        new_parent_state      = parent_state;
                        goto Found;
                    }
                }
            }
            
            // now that we've chosen the new states, we need to assign them
            Found:
            
                static_cast<CharacterEventDiscrete*>(nodeParentState[site_index])->setState(new_parent_state);
                static_cast<CharacterEventDiscrete*>(brotherParentState[site_index])->setState(new_parent_state);
                static_cast<CharacterEventDiscrete*>(parentChildState[site_index])->setState(new_parent_state);
            
                static_cast<CharacterEventDiscrete*>(uncleParentState[site_index])->setState(new_grandparent_state);
                static_cast<CharacterEventDiscrete*>(parentParentState[site_index])->setState(new_grandparent_state);
                static_cast<CharacterEventDiscrete*>(grandparentChildState[site_index])->setState(new_grandparent_state);
            
        }
        
    }
    else
    {
     
        // update the transition probs
        std::vector<double> rf = c->getRootFrequencies();
        
        // for each site...
        for (size_t site_index = 0; site_index < num_sites; ++site_index)
        {
            
            size_t node_state              = static_cast<CharacterEventDiscrete*>(nodeChildState[site_index])->getState();
            size_t brother_state           = static_cast<CharacterEventDiscrete*>(brotherChildState[site_index])->getState();
            size_t uncle_state             = static_cast<CharacterEventDiscrete*>(uncleChildState[site_index])->getState();
            
            double sum = 0.0;
            
            std::vector<double> state_probs(num_states * num_states);
            
            // for each grandparent state...
            for (size_t grandparent_state = 0; grandparent_state < num_states; ++grandparent_state)
            {
                
                // for each parent state...
                for (size_t parent_state = 0; parent_state < num_states; ++parent_state)
                {
                    double p = parent_tp_matrix[grandparent_state][parent_state] * node_tp_matrix[parent_state][node_state] * brother_tp_matrix[parent_state][brother_state] * uncle_tp_matrix[grandparent_state][uncle_state] * rf[grandparent_state];
                    sum += p;
                    state_probs[grandparent_state * num_states + parent_state] = p;
                }
                
            }
            
            // generate a random uniform number
            double u = GLOBAL_RNG->uniform01();
            
            size_t new_grandparent_state = 0;
            size_t new_parent_state      = 0;
            
            // sample a new state
            for (size_t grandparent_state = 0; grandparent_state < num_states; ++grandparent_state)
            {
                for (size_t parent_state = 0; parent_state < num_states; ++parent_state)
                {
                    u -= (state_probs[grandparent_state * num_states + parent_state] / sum);
                    if ( u < 0.0 )
                    {
                        new_grandparent_state = grandparent_state;
                        new_parent_state      = parent_state;
                        goto FoundRoot;
                    }
                }
            }
            
            // now that we've chosen the new states, we need to assign them
            FoundRoot:
            
                static_cast<CharacterEventDiscrete*>(nodeParentState[site_index])->setState(new_parent_state);
                static_cast<CharacterEventDiscrete*>(brotherParentState[site_index])->setState(new_parent_state);
                static_cast<CharacterEventDiscrete*>(parentChildState[site_index])->setState(new_parent_state);
            
                static_cast<CharacterEventDiscrete*>(uncleParentState[site_index])->setState(new_grandparent_state);
                static_cast<CharacterEventDiscrete*>(parentParentState[site_index])->setState(new_grandparent_state);
                static_cast<CharacterEventDiscrete*>(grandparentChildState[site_index])->setState(new_grandparent_state);

        }

    }

}

template<class charType>
void RevBayesCore::NarrowExchangeCharacterHistoryProposal<charType>::sampleNodeCharacters( TopologyNode* node )
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
        rm.calculateTransitionProbabilities(node_age, left_age,  left_rate,  brother_tp_matrix);
        rm.calculateTransitionProbabilities(node_age, right_age, right_rate, uncle_tp_matrix);
        
        // states for conditional sampling probs
        const std::vector<CharacterEvent*>& leftChildState  = histories[left_child.getIndex()]->getChildCharacters();
        const std::vector<CharacterEvent*>& rightChildState = histories[right_child.getIndex()]->getChildCharacters();
        
        // states to update
        std::vector<CharacterEvent*>& nodeChildState   = histories[node->getIndex()]->getChildCharacters();
        std::vector<CharacterEvent*>& leftParentState  = histories[left_child.getIndex()]->getParentCharacters();
        std::vector<CharacterEvent*>& rightParentState = histories[right_child.getIndex()]->getParentCharacters();
        
        if ( node->isRoot() == false )
        {
            double parent_age   = node->getParent().getAge();
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
                    double p = node_tp_matrix[ancS][i] * brother_tp_matrix[i][desS1] * uncle_tp_matrix[i][desS2];
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
                    double p = rf[i] * brother_tp_matrix[i][desS1] * uncle_tp_matrix[i][desS2];
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
void RevBayesCore::NarrowExchangeCharacterHistoryProposal<charType>::setRateGenerator(const TypedDagNode<RateGenerator> *d)
{
    
    rate_generator = d;
    num_states = rate_generator->getValue().getNumberOfStates();
    
    node_proposal->setRateGenerator( rate_generator );
    parent_proposal->setRateGenerator( rate_generator );
    uncle_proposal->setRateGenerator( rate_generator );
    brother_proposal->setRateGenerator( rate_generator );
    grandparent_proposal->setRateGenerator( rate_generator );
    
    node_tp_matrix        = TransitionProbabilityMatrix( num_states );
    brother_tp_matrix     = TransitionProbabilityMatrix( num_states );
    uncle_tp_matrix       = TransitionProbabilityMatrix( num_states );
    parent_tp_matrix      = TransitionProbabilityMatrix( num_states );
    grandparent_tp_matrix = TransitionProbabilityMatrix( num_states );
 
}


template<class charType>
void RevBayesCore::NarrowExchangeCharacterHistoryProposal<charType>::setRateGenerator(const TypedDagNode<RateGeneratorSequence> *d)
{
    
    rate_generator_sequence = d;
    num_states = rate_generator_sequence->getValue().getNumberOfStates();
    
    node_proposal->setRateGenerator( rate_generator_sequence );
    parent_proposal->setRateGenerator( rate_generator_sequence );
    uncle_proposal->setRateGenerator( rate_generator_sequence );
    brother_proposal->setRateGenerator( rate_generator_sequence );
    grandparent_proposal->setRateGenerator( rate_generator_sequence );

    node_tp_matrix        = TransitionProbabilityMatrix( num_states );
    brother_tp_matrix     = TransitionProbabilityMatrix( num_states );
    uncle_tp_matrix       = TransitionProbabilityMatrix( num_states );
    parent_tp_matrix      = TransitionProbabilityMatrix( num_states );
    grandparent_tp_matrix = TransitionProbabilityMatrix( num_states );
    
}

/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
template<class charType>
void RevBayesCore::NarrowExchangeCharacterHistoryProposal<charType>::swapNodeInternal(DagNode *oldN, DagNode *newN)
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
    
    node_proposal->swapNodeInternal(oldN, newN);
    parent_proposal->swapNodeInternal(oldN, newN);
    uncle_proposal->swapNodeInternal(oldN, newN);
    brother_proposal->swapNodeInternal(oldN, newN);
    grandparent_proposal->swapNodeInternal(oldN, newN);
    
}


template<class charType>
void RevBayesCore::NarrowExchangeCharacterHistoryProposal<charType>::setProposalTuningParameter(double tp)
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
void RevBayesCore::NarrowExchangeCharacterHistoryProposal<charType>::tune( double rate )
{
    
    // nothing to tune
    
}

#endif
