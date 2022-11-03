#ifndef NodeTimeSlideUniformCharacterHistoryProposal_H
#define NodeTimeSlideUniformCharacterHistoryProposal_H

#include "BranchHistoryDiscrete.h"
#include "CharacterEventDiscrete.h"
#include "DistributionUniform.h"
#include "PathRejectionSampleProposal.h"
#include "Proposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RateGenerator.h"
#include "RateGeneratorSequence.h"
#include "RbConstants.h"
#include "RbException.h"
#include "StochasticNode.h"
#include "Tree.h"
#include "TreeHistoryCtmc.h"
#include "TreeUtilities.h"
#include "TypedDagNode.h"

#include <cmath>
#include <iostream>
#include <string>

namespace RevBayesCore {
    
    /**
     * The node-age slide proposal operator using a Uniform distribution.
     *
     * This node-age proposal is a Uniform-sliding proposal on rooted subtrees without changing the topology.
     * That is, we pick a random node which is not the root.
     * Then, we pick an age between the parent and the oldest sampled descendant drawn from a Uniform distribution centered around the current age.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2012-07-12, version 1.0
     *
     */
    
    template<class charType>
    class NodeTimeSlideUniformCharacterHistoryProposal : public Proposal {
        
    public:
        NodeTimeSlideUniformCharacterHistoryProposal( StochasticNode<Tree> *t, StochasticNode<AbstractHomologousDiscreteCharacterData>* c);                                               //!<  constructor
        
        // Basic utility functions
        void                                                        cleanProposal(void);                                        //!< Clean up proposal
        NodeTimeSlideUniformCharacterHistoryProposal*               clone(void) const;                                          //!< Clone object
        double                                                      doProposal(void);                                           //!< Perform proposal
        const std::string&                                          getProposalName(void) const;                                //!< Get the name of the proposal for summary printing
        void                                                        prepareProposal(void);                                      //!< Prepare the proposal
        void                                                        printParameterSummary(std::ostream &o, bool name_only) const;               //!< Print the parameter summary
        void                                                        sampleNodeCharacters(TopologyNode* node);                       //!< Sample the characters at the node
        void                                                        setRateGenerator(const TypedDagNode<RateGenerator> *d);
        void                                                        setRateGenerator(const TypedDagNode<RateGeneratorSequence> *d);
        void                                                        undoProposal(void);                                         //!< Reject the proposal
        
    protected:
        
        void                                                        swapNodeInternal(DagNode *oldN, DagNode *newN);             //!< Swap the DAG nodes on which the Proposal is working on
        
        
    private:
        
        // parameters
        StochasticNode<Tree>*                                       tree;                                                           //!< The tree the Proposal is working on
        StochasticNode<AbstractHomologousDiscreteCharacterData>*    ctmc;                                                           //!< The character history the Proposal is working on
        const TypedDagNode<RateGenerator>*                          rate_generator;
        const TypedDagNode<RateGeneratorSequence>*                  rate_generator_sequence;

        // dimensions
        size_t                                                      num_states;

        // proposals
//        NodeRejectionSampleProposal<charType>*                      node_proposal;
        PathRejectionSampleProposal<charType>*                      node_proposal;
        PathRejectionSampleProposal<charType>*                      left_proposal;
        PathRejectionSampleProposal<charType>*                      right_proposal;
        
        TransitionProbabilityMatrix                                 node_tp_matrix;
        TransitionProbabilityMatrix                                 left_tp_matrix;
        TransitionProbabilityMatrix                                 right_tp_matrix;

        // stored objects to undo proposal
        TopologyNode*                                               stored_node;
        double                                                      stored_age;
        std::vector<size_t>                                         stored_node_states;
        
    };
    
}

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
template<class charType>
RevBayesCore::NodeTimeSlideUniformCharacterHistoryProposal<charType>::NodeTimeSlideUniformCharacterHistoryProposal( StochasticNode<Tree> *t, StochasticNode<AbstractHomologousDiscreteCharacterData>* c ) : Proposal(),
    tree( t ),
    ctmc( c ),
    rate_generator( NULL ),
    rate_generator_sequence( NULL ),
    num_states( 2 ),
    node_tp_matrix( 2 ),
    left_tp_matrix( 2 ),
    right_tp_matrix( 2 )
{
    
    // tell the base class to add the node
    addNode( tree );
    addNode( ctmc );

    // make the proposal
    node_proposal  = new PathRejectionSampleProposal<charType>( c );
    left_proposal  = new PathRejectionSampleProposal<charType>( c );
    right_proposal = new PathRejectionSampleProposal<charType>( c );
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
template<class charType>
void RevBayesCore::NodeTimeSlideUniformCharacterHistoryProposal<charType>::cleanProposal( void )
{
    
    node_proposal->cleanProposal();
    left_proposal->cleanProposal();
    right_proposal->cleanProposal();
    
    stored_node_states.clear();
    
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
template<class charType>
RevBayesCore::NodeTimeSlideUniformCharacterHistoryProposal<charType>* RevBayesCore::NodeTimeSlideUniformCharacterHistoryProposal<charType>::clone( void ) const
{
    
    return new NodeTimeSlideUniformCharacterHistoryProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
template<class charType>
const std::string& RevBayesCore::NodeTimeSlideUniformCharacterHistoryProposal<charType>::getProposalName( void ) const
{
    static std::string name = "NodeTimeSlideUniformCharacterHistory";
    
    return name;
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
template<class charType>
double RevBayesCore::NodeTimeSlideUniformCharacterHistoryProposal<charType>::doProposal( void )
{
    
    // get random number generator
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    // retrieve the tree
    Tree& tau = tree->getValue();
    
    // pick a random node which is neither the root nor the direct descendant of the root
    TopologyNode* node;
    do {
        double u = rng->uniform01();
        size_t index = size_t( std::floor(tau.getNumberOfNodes() * u) );
        node = &tau.getNode(index);
    } while ( node->isRoot() || node->isTip() );
    
    TopologyNode& parent = node->getParent();
    
    // we need to work with the times
    double parent_age = parent.getAge();
    double my_age     = node->getAge();
    double child_age  = node->getChild( 0 ).getAge();
    if ( child_age < node->getChild( 1 ).getAge())
    {
        child_age = node->getChild( 1 ).getAge();
    }

    // retrieve the character history
    TreeHistoryCtmc<charType>* p = dynamic_cast< TreeHistoryCtmc<charType>* >( &ctmc->getDistribution() );
    if ( p == NULL )
    {
        throw RbException("Failed cast.");
    }
    size_t num_sites = p->getNumberOfSites();
    
    // now we store all necessary values
    stored_node = node;
    stored_age  = my_age;
    
    stored_node_states.clear();
    const std::vector<CharacterEvent*>& node_states = p->getHistory( *node ).getChildCharacters();
    stored_node_states.resize(num_sites, 0);
    for (size_t site_index = 0; site_index < num_sites; ++site_index)
    {
        size_t s = static_cast<CharacterEventDiscrete*>(node_states[site_index])->getState();
        stored_node_states[site_index] = s;
    }

    // assign the nodes for the path proposals
    node_proposal->assignNode( node );
    left_proposal->assignNode( &node->getChild( 0 ) );
    right_proposal->assignNode( &node->getChild( 1 ) );
    
    // prepare the path proposals
    node_proposal->prepareProposal();
    left_proposal->prepareProposal();
    right_proposal->prepareProposal();

    // draw new age for the chosen node
    double my_new_age = (parent_age - child_age) * rng->uniform01() + child_age;
    
    // set the age
    tau.getNode(node->getIndex()).setAge( my_new_age );
    
    // sample new node states for the node
    sampleNodeCharacters( node );

    // sample new branch histories
    double ln_proposal_probability = 0.0;
    ln_proposal_probability += node_proposal->doProposal();
    ln_proposal_probability += left_proposal->doProposal();
    ln_proposal_probability += right_proposal->doProposal();

    return ln_proposal_probability;
    
}


/**
 *
 */
template<class charType>
void RevBayesCore::NodeTimeSlideUniformCharacterHistoryProposal<charType>::prepareProposal( void )
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
void RevBayesCore::NodeTimeSlideUniformCharacterHistoryProposal<charType>::printParameterSummary(std::ostream &o, bool name_only) const
{
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
template<class charType>
void RevBayesCore::NodeTimeSlideUniformCharacterHistoryProposal<charType>::undoProposal( void )
{

    // retrieve the tree
    Tree& tau = tree->getValue();

    // reset the node age
    tau.getNode(stored_node->getIndex()).setAge( stored_age );
//    stored_node->setAge( stored_age );

    // retrieve the history
    TreeHistoryCtmc<charType>* p = dynamic_cast< TreeHistoryCtmc<charType>* >(&ctmc->getDistribution());
    if ( p == NULL )
    {
        throw RbException("Failed cast.");
    }
    size_t num_sites = p->getNumberOfSites();
    const std::vector<BranchHistory*>& histories = p->getHistories();

    // restore node state
    std::vector<CharacterEvent*>& node_child_state   = histories[stored_node->getIndex()]->getChildCharacters();
    std::vector<CharacterEvent*>& left_parent_state  = histories[stored_node->getChild(0).getIndex()]->getParentCharacters();
    std::vector<CharacterEvent*>& right_parent_state = histories[stored_node->getChild(1).getIndex()]->getParentCharacters();
    for (size_t site_index = 0; site_index < num_sites; ++site_index)
    {
        size_t s = stored_node_states[site_index];
        static_cast<CharacterEventDiscrete*>(node_child_state[site_index])->setState(s);
        static_cast<CharacterEventDiscrete*>(left_parent_state[site_index])->setState(s);
        static_cast<CharacterEventDiscrete*>(right_parent_state[site_index])->setState(s);
    }

    // restore path state
    node_proposal->undoProposal();
    left_proposal->undoProposal();
    right_proposal->undoProposal();

}

template<class charType>
void RevBayesCore::NodeTimeSlideUniformCharacterHistoryProposal<charType>::sampleNodeCharacters( TopologyNode* node )
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
void RevBayesCore::NodeTimeSlideUniformCharacterHistoryProposal<charType>::setRateGenerator(const TypedDagNode<RateGenerator> *d)
{
    
    rate_generator = d;
    num_states = rate_generator->getValue().getNumberOfStates();
    
    node_proposal->setRateGenerator( rate_generator );
    left_proposal->setRateGenerator( rate_generator );
    right_proposal->setRateGenerator( rate_generator );
    
    node_tp_matrix  = TransitionProbabilityMatrix( num_states );
    left_tp_matrix  = TransitionProbabilityMatrix( num_states );
    right_tp_matrix = TransitionProbabilityMatrix( num_states );
    
}


template<class charType>
void RevBayesCore::NodeTimeSlideUniformCharacterHistoryProposal<charType>::setRateGenerator(const TypedDagNode<RateGeneratorSequence> *d)
{
    
    rate_generator_sequence = d;
    num_states = rate_generator_sequence->getValue().getNumberOfStates();
    
    node_proposal->setRateGenerator( rate_generator_sequence );
    left_proposal->setRateGenerator( rate_generator_sequence );
    right_proposal->setRateGenerator( rate_generator_sequence );
    
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
void RevBayesCore::NodeTimeSlideUniformCharacterHistoryProposal<charType>::swapNodeInternal(DagNode *oldN, DagNode *newN)
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
    left_proposal->swapNodeInternal(oldN, newN);
    right_proposal->swapNodeInternal(oldN, newN);
    
}


#endif
